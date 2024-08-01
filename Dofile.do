cd "C:/Users/Usuario/OneDrive/Escritorio/3ro GANE/Mercados/Stata/1"
use "data.dta", clear

*INTRODUCCIÓN:
describe
**Seleccionamos variables relevantes
keep gasto_azucar renta_pc gasto_cafenormal gasto_cafecapsulas gasto_lecheentera gasto_lechedescremada sexosp paisnacsp estudiossp

********************************************************************************
********************************************************************************

*PARTE 1: ANÁLISIS ESTADÍSTICO:
********************************************************************************
** Variables principales

*** Gasto en azucar del sustentador principal de las familias consumidoras
sum gasto_azucar, detail
****Corregimos valores outliers
replace gasto_azucar=. if gasto_azucar<0.1
egen meanazu=mean(gasto_azucar)
egen stdazu=sd(gasto_azucar)
gen zazu=meanazu+(2.5*stdazu)
keep if gasto_azucar<zazu
****Visualizamos resultados
sum gasto_azucar, detail
histogram gasto_azucar

*** Renta
gen renta=(renta_pc/12)/100
sum renta, detail
****Corregimos valores outliers
replace renta=. if renta<1.19
egen meanren=mean(renta)
egen stdren=sd(renta)
gen zren=meanren+(2.5*stdren)
keep if renta<zren
****Observamos distribución
sum renta, detail
histogram renta

*** Introducimos Relación inicial
graph twoway (scatter gasto_azucar renta) (lpoly gasto_azucar renta) (lowess gasto_azucar renta), title("Relación Gasto Azúcar y Renta") saving(1.gph,replace)

***Definimos globales
global y gasto_azucar
global x1 renta

********************************************************************************
** Variables complementarias
***Gasto en cafe y leche
gen caf=gasto_cafenormal+gasto_cafecapsulas+gasto_lecheentera+gasto_lechedescremada
egen c33 = pctile(caf), p(33)
egen c66 = pctile(caf), p(66)

gen gcyl=1 if caf<c33
replace gcyl=2 if caf>c33 & caf<c66
replace gcyl=3 if caf>c66
label define cf 1 "Menor Consumo" 2 "Consumo Intermedio" 3 "Mayor Consumo"
label value gcyl cf
corr gasto_azucar $y
reg gasto_azucar renta caf, robust

***definimos globales
global x2 i.gcyl

********************************************************************************
** Variables por individuos
***Educación
gen educ=1 if estudiossp==1 | estudiossp==2
replace educ=2 if estudiossp==3 | estudiossp==4
replace educ=3 if estudiossp==5
replace educ=4 if estudiossp==6 | estudiossp==7 | estudiossp==8

label define ed 1 "Educacion primaria o inferior" 2 "Educacion secundaria, Bachiller o FP" 3 "Grado superior" 4 "Graduado universitario o superior"
label value educ ed
reg gasto_azucar renta i.gcyl i.educ, robust

***definimos globales
global x3  i.educ

********************************************************************************
********************************************************************************

*PARTE 2: MODELIZACIÓN:
**Estimación formas funcionales
*** Lin
qui reg $y $x1 $x2 $x3, robust
estimate store mlin

*** Lin-Log
gen lrenta=log(renta)
qui reg $y lrenta $x2 $x3, robust
estimate store mlinlog

*** Log-Lin
gen lgasazu=log(gasto_azucar)
qui reg lgasazu $x1 $x2 $x3, robust
estimate store mloglin

*** Log-Log
qui reg lgasazu lrenta $x2 $x3, robust
estimate store mloglog

*** Sq
qui reg $y $x1 c.renta#c.renta $x2 $x3, robust
estimate store msq

*** Representamos resultados
estimate table mlin mlinlog mloglin mloglog msq, stats(N r2_a)

********************************************************************************
**Obtención efectos marginales
***Sobre la media
**** Lin
qui reg $y $x1 $x2 $x3, robust
scalar blin=e(b)[1,1]
gen EMGlin=blin
**** Lin-Log
qui reg $y lrenta $x2 $x3, robust
scalar blinlog=e(b)[1,1]
egen meanlogren=mean(lrenta)
gen EMGlinlog=blinlog/meanlogren
**** Log-Lin
qui reg lgasazu $x1 $x2 $x3, robust
scalar bloglin=e(b)[1,1]
egen meanlogazu=mean(lgasazu)
gen EMGloglin=bloglin*meanlogazu
**** Log-Log
qui reg lgasazu lrenta $x2 $x3, robust
scalar bloglog=e(b)[1,1]
gen EMGloglog=bloglog*(meanlogazu/meanlogren)
**** Sq
qui reg $y $x1 c.renta#c.renta $x2 $x3, robust
scalar b1sq=e(b)[1,1]
scalar b2sq=e(b)[1,2]
egen meanx=mean(renta)
gen EMGsq=b1sq+(2*b2sq*meanx)
**** Representamos resultados
tabstat EMGlin EMGlinlog EMGloglin EMGloglog EMGsq
//Para ver valores de scalar: di

***Sobre p10 
**** Lin-Log
egen lr10 = pctile(lrenta), p(10)
gen EMG10linlog=blinlog/lr10
**** Log-Lin
egen la10=pctile(lgasazu), p(10)
gen EMG10loglin=bloglin*la10
**** Log-Log
gen EMG10loglog=bloglog*(la10/lr10)
**** Sq
egen r10=pctile(renta), p(10)
gen EMG10sq=b1sq+(2*b2sq*r10)
**** Representamos resultados
tabstat EMGlin EMG10linlog EMG10loglin EMG10loglog EMG10sq

***Sobre p50 (Mediana)
**** Lin-Log
egen lr50 = pctile(lrenta), p(50)
gen EMG50linlog=blinlog/lr50
**** Log-Lin
egen la50=pctile(lgasazu), p(50)
gen EMG50loglin=bloglin*la50
**** Log-Log
gen EMG50loglog=bloglog*(la50/lr50)
**** Sq
egen r50=pctile(renta), p(50)
gen EMG50sq=b1sq+(2*b2sq*r50)
**** Representamos resultados
tabstat EMGlin EMG50linlog EMG50loglin EMG50loglog EMG50sq


***Sobre p90
**** Lin-Log
egen lr90 = pctile(lrenta), p(90)
gen EMG90linlog=blinlog/lr90
**** Log-Lin
egen la90=pctile(lgasazu), p(90)
gen EMG90loglin=bloglin*la90
**** Log-Log
gen EMG90loglog=bloglog*(la90/lr90)
**** Sq
egen r90=pctile(renta), p(90)
gen EMG90sq=b1sq+(2*b2sq*r90)
**** Representamos resultados
tabstat EMGlin EMG90linlog EMG90loglin EMG90loglog EMG90sq


********************************************************************************
**Obtención elasticidades
***Sobre la media
**** Lin
egen meany=mean(gasto_azucar)
gen Elin=blin*(meanx/meany)
**** Lin-Log
gen Elinlog=blinlog/meany
**** Log-Lin
gen Eloglin=bloglin*meanx
**** Log-Log
gen Eloglog=bloglog
**** Sq
gen Esq=EMGsq*(meanx/meany)
**** Representamos resultados
tabstat Elin Elinlog Eloglin Eloglog Esq

***Sobre p10 
**** Lin-Log
egen y10=pctile(gasto_azucar), p(10)
egen x10=pctile(renta), p(10)
gen E10lin=blin*(x10/y10)
**** Lin-Log
gen E10linlog=blinlog/y10
**** Log-Lin
gen E10loglin=bloglin*x10
**** Sq
gen E10sq=EMGsq*(x10/y10)
**** Representamos resultados
tabstat E10lin E10linlog E10loglin Eloglog E10sq

***Sobre p50 (Mediana)
**** Lin-Log
egen y50=pctile(gasto_azucar), p(50)
egen x50=pctile(renta), p(50)
gen E50lin=blin*(x50/y50)
**** Lin-Log
gen E50linlog=blinlog/y50
**** Log-Lin
gen E50loglin=bloglin*x50
**** Sq
gen E50sq=EMGsq*(x50/y50)
**** Representamos resultados
tabstat E50lin E50linlog E50loglin Eloglog E50sq

***Sobre p90
**** Lin-Log
egen y90=pctile(gasto_azucar), p(90)
egen x90=pctile(renta), p(90)
gen E90lin=blin*(x90/y90)
**** Lin-Log
gen E90linlog=blinlog/y90
**** Log-Lin
gen E90loglin=bloglin*x90
**** Sq
gen E90sq=EMGsq*(x90/y90)
**** Representamos resultados
tabstat E90lin E90linlog E90loglin Eloglog E90sq

********************************************************************************
**Selección forma funcional

***Sobre la media
**** Lin
qui reg $y $x1 $x2 $x3, robust
predict plin
*** Lin-Log
qui reg $y lrenta $x2 $x3, robust
predict plinlog
**** Sq
qui reg $y $x1 c.renta#c.renta $x2 $x3, robust
predict psq
**** Log-Lin
qui reg lgasazu $x1 $x2 $x3, robust
predict loglin
predict resloglin, res
gen expresloglin=exp(resloglin)
egen fcloglin=mean(expresloglin)
gen ploglin=exp(loglin)*fcloglin
**** Log-Log
qui reg lgasazu lrenta $x2 $x3, robust
predict loglog
predict resloglog, res
gen expresloglog=exp(resloglog)
egen fcloglog=mean(expresloglog)
gen ploglog=exp(loglog)*fcloglog
**** Comprobamos resultados
tabstat meany plin plinlog ploglin ploglog psq
**** Diferencias con respecto a la media
gen diflin=plin-meany
gen diflinlog=plinlog-meany
gen difloglin=ploglin-meany
gen difloglog=ploglog-meany
gen difsq=psq-meany
tabstat diflin diflinlog difloglin difloglog difsq

***Sobre p25
**** Hallamos valores
egen p25y=pctile(gasto_azucar), p(25)
egen p25lin=pctile(plin), p(25)
egen p25linlog=pctile(plinlog), p(25)
egen p25loglin=pctile(ploglin), p(25)
egen p25loglog=pctile(ploglog), p(25)
egen p25sq=pctile(psq), p(25)
**** Comprobamos resultados
tabstat p25y p25lin p25linlog p25loglin p25loglog p25sq
**** Diferencias con respecto a p25
gen dif25lin=p25lin-p25y
gen dif25linlog=p25linlog-p25y
gen dif25loglin=p25loglin-p25y
gen dif25loglog=p25loglog-p25y
gen dif25sq=p25sq-p25y
tabstat dif25lin dif25linlog dif25loglin dif25loglog dif25sq

***Sobre p50 (Mediana)
**** Hallamos valores
egen p50y=pctile(gasto_azucar), p(50)
egen p50lin=pctile(plin), p(50)
egen p50linlog=pctile(plinlog), p(50)
egen p50loglin=pctile(ploglin), p(50)
egen p50loglog=pctile(ploglog), p(50)
egen p50sq=pctile(psq), p(50)
**** Comprobamos resultados
tabstat p50y p50lin p50linlog p50loglin p50loglog p50sq
**** Diferencias con respecto a p50
gen dif50lin=p50lin-p50y
gen dif50linlog=p50linlog-p50y
gen dif50loglin=p50loglin-p50y
gen dif50loglog=p50loglog-p50y
gen dif50sq=p50sq-p50y
tabstat dif50lin dif50linlog dif50loglin dif50loglog dif50sq

***Sobre p75
**** Hallamos valores
egen p75y=pctile(gasto_azucar), p(75)
egen p75lin=pctile(plin), p(75)
egen p75linlog=pctile(plinlog), p(75)
egen p75loglin=pctile(ploglin), p(75)
egen p75loglog=pctile(ploglog), p(75)
egen p75sq=pctile(psq), p(75)
**** Comprobamos resultados
tabstat p75y p75lin p75linlog p75loglin p75loglog p75sq
**** Diferencias con respecto a p25
gen dif75lin=p75lin-p75y
gen dif75linlog=p75linlog-p75y
gen dif75loglin=p75loglin-p75y
gen dif75loglog=p75loglog-p75y
gen dif75sq=p75sq-p75y
tabstat dif75lin dif75linlog dif75loglin dif75loglog dif75sq

********************************************************************************
********************************************************************************

*PARTE 3: ENFOQUE ALTERNATIVO
** Curva Engel
***Redefinimos renta
scalar  tp=-b1sq/(2*b2sq)
gen x11=renta
replace x11=0 if renta>tp
gen x12=renta
replace x12=0 if renta<tp
reg $y x11 x12 $x2 $x3, robust
estimate store engel
test x11=x12

****Representamos coeficientes
coefplot(engel), keep(x11 x12) yline(0) vertical title("Efectos Marginales Renta-Consumo de Azúcar") saving(engel.gph, replace)

********************************************************************************
** Segmentos por características
***Concretamos la variable renta
gen edprisec=0
replace edprisec=1 if educ==1 | educ==2
gen renta_edprisec=edprisec*renta

gen gradsup=0
replace gradsup=1 if educ==3
gen renta_gradsup=gradsup*renta

gen uni=0
replace uni=1 if educ==4
gen renta_uni=uni*renta

global x32 gradsup uni

***Estimamos
reg $y renta_edprisec renta_gradsup renta_uni $x2 $x32, robust
estimate store seg
test  renta_edprisec -renta_gradsup=0
test  renta_gradsup-renta_uni=0, accumulate

****Representamos coeficientes
coefplot(seg), keep(renta_edprisec renta_gradsup renta_uni) yline(0) vertical title("Efectos Marginales Renta-Consumo de Azúcar") saving(seg.gph, replace)