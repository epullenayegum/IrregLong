# IrregLong
R package for analysing longitudinal data subject to irregular observation

[![Travis build status](https://travis-ci.org/epullenayegum/IrregLong.svg?branch=master)](https://travis-ci.org/epullenayegum/IrregLong)

Why should I use it?
Irregular follow-up in longitudinal data can lead to biased inferences if the visit process is related to the outcome of interest. For example, in a study of HIV positive mothers in Kenya, women visited a clinic according to a pre-defined schedule, but also visited as needed (for example, due to illness). The estimated prevalence of pneumonia in the first year postpartum was 2.89% when the visit process was ignored, but almost halved (estimated prevalence 1.48%) when the visit process was accounted for.

How do I use it?
The main command in the package is iiwgee; this fits an inverse-intensity weighted GEE and will result in consistent estimates of regression coefficients provided the visit and outcome processes are conditionally independent given a set of observed covariates.

iiwgee consists of two models; one for the outcome of interest, and the other for the visit process. For example, in a study of juvenile dermatomyositis, the outcome of interest was a disease activity score (DAS) and the visit process was described through the following model for the visit intensity $\lambda$ at time $t$:
\[\lambda_i(t)=\lambda_0(t)exp(\gamma log(DASlag_i(t)+1)) \]
where $DASlag_i(t)$ is the last observed disease activiy score for patient $i$ at time $t$. Suppose we are interested in modelling the mean disease activity as a linear function of time. The syntax for fitting the inverse intensity weighted GEE is then

```{r,eval=FALSE}
iiwgee(DAS~time,Surv(time.lag,time,event)~log(DAS.lag+1),data=data,lagvars=("time","DAS"),id="id",time="time",event="event",invariant="id",maxfu=4,first=TRUE)
```
Please see the package documentation for more details.

How do I get it?
Install the devtools package from CRAN
Load the package (type: library(devtools))
Type: install_github("epullenayegum/IrregLong/IrregLong")
