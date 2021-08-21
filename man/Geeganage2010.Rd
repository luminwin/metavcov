\name{Geeganage2010}
\alias{Geeganage2010}
\docType{data}
\title{17 studies of multivariate effect sizes reported by Geeganage et al. (2010) }
\description{
This dataset includes 17 studies of multivariate effect sizes with four different outcomes reported by Geeganage and Bath (2010).}
  \usage{data(Geeganage2010)}
\details{
In a meta-analysis, Geeganage and Bath (2010) studied whether blood pressure (BP) should be actively altered during the acute phase of stroke, and assessed the effect of multiple vasoactive drugs on BP in acute stroke. Selection criteria included: Randomized trials of interventions that would be expected, on pharmacological grounds, to alter BP in patients within one week of the onset of acute stroke. There were four outcomes: systolic blood pressure (SBP, in mHg), diastolic blood pressure (DBP, in mHg), death (D), and death or disability (DD).

\tabular{ll}{
ID: \tab ID for each study included\cr
ft_D \tab	Number of early death within 1 month (D) in "1 Drug" Group\cr
fc_D \tab	Number of D in "control " Group \cr
nt_D \tab	Number of people in "1 Drug" Group reporting D status \cr
nc_D \tab	Number of people in "control " Group reporting D status \cr
OR_D \tab	Odds Ratio of D for "1 Drug" versus "control" group \cr
ft_DD \tab	Number of early death or deterioration within 1 month (DD) in "1 Drug" Group \cr
fc_DD \tab	Number of early DD in "control " Group \cr
nt_DD \tab	Number of people in "1 Drug" Group reporting DD status \cr
nc_DD \tab	Number of people in "control " Group reporting DD status \cr
OR_DD \tab	Odds Ratio of DD for "1 Drug" versus "control" group \cr
nt_SBP \tab	Number of people in "1 Drug" Group reporting Systolic blood pressure (SBP) status \cr
nc_SBP \tab	Number of people in "control " Group reporting SBP status \cr
MD_SBP \tab	Mean Difference of SBP for "1 Drug" versus "control" group \cr
sdt_SBP \tab	Standard Deviation of SBP in "1 Drug" Group \cr
sdc_SBP \tab	Standard Deviation of SBP in "control " Group \cr
nt_DBP \tab	Number of people in "1 Drug" Group reporting Diastolic blood pressure (DBP) status \cr
nc_DBP \tab	Number of people in "control " Group reporting DBP status \cr
MD_DBP \tab	Mean Difference of DBP for "1 Drug" versus "control" group \cr
sdt_DBP \tab	Standard Deviation of DBP in "1 Drug" Group \cr
sdc_DBP \tab	Standard Deviation of DBP in "control " Group \cr
SMD_SBP \tab	Standardized Mean Difference of SBP for "1 Drug" versus "control" group \cr
SMD_DBP \tab	Standardized Mean Difference of DBP for "1 Drug" versus "control" group \cr}
}
\examples{
  data(Geeganage2010)
}
\source{
Geeganage, C., & Bath, P. M. (2010). Vasoactive drugs for acute stroke. \emph{The Cochrane Library}.
}


