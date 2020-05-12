[TOC]

# STAT 507 PennState

 https://online.stat.psu.edu/stat507/ 



### 1.4 - Hypotheses in Epidemiology, Designs and Populations

Types of studies

- Case Study (describing one person with the condition, a case)
- Case Series (series of cases)
- Ecological Study (analysis of group statistics..for example, comparing rates of disease between two countries)
- Cross-Sectional Study (assessing individuals at one time, such as a survey)
- Case-Control Study (studying those with the condition vs. those without)
- Cohort Study (following subjects over time to study the initiation and progression of a condition)

### 2.2 - Measures of Disease Frequency

The commonly used measures of *incidence* and *prevalence* can be distinguished by differences in the time of disease onset 

* **Incidence** is a count of *new cases* of the disease (or outcome). 
* **Prevalence**, on the other hand, counts *both new and existing cases* of the disease. 

Epidemiologic measures of disease frequency are of 5 types:

1. Count
2. Proportion  A/(A+B) 
3. Ratio A/B:  A ratio *as a measure of disease frequency* is used infrequently, in special situations 
4. Rate:  An epidemiologic rate will contain the following: disease frequency (numerator), unit of population size, and the time period during which the event occurred 
5. Risk:   the probability of an individual meeting the case definition (person-time rate). Risk is dependent upon time. 

### 2.3 - Disease Occurrence

**Incidence**: counts *new* cases of the disease (or outcome)

* A summary incidence rate can estimate risk (e.g., probability of disease in an individual) if risk is constant across the summarized groups. 
* It is often expressed as a proportion of those at risk. The denominator includes all persons at-risk for the disease or condition, i.e. disease-free or condition-free individuals in the population at the start of the time period 
* Incidence can also be expressed in terms of person-time at risk. 
* Rates are usually expressed per 100, 1,000, or 100,000 persons. 
* Two Common Measures of Incidence:
  * cumulative incidence rate  =  incidence rate :  the number of persons who newly experience the disease or studied outcome during a specified period of time divided by the (average) total population at risk 
  *  Incidence density rate =  incidence rate; person-time rate, force of morbidity/mortality, hazard rate,disease intensity : the number of persons who newly experience the outcome during a specified period of time divided by the sum of the time that each member of the population is at-risk 

**Prevalence**: counts *new and existing* cases of the disease (or outcome)

* Diseases with a long duration will be more prevalent than those with shorter duration ：Prevalence = Incidence × Average Duration 
* Prevalence is often expressed after multiplication by 100 (%), 1000 or 100,000. 
* Prevalence is a proportion, usually reflecting the proportion with a disease at a particular time.  
* The **prevalence pool** is the subset of the population with the condition of interest.  
* **Point prevalence**: prevalence of condition of interest at a specific time. 
* **Period prevalence**: prevalence of outcome of interest during a specified period of time.  

### 3.1 - Exposure Concepts

**Exposure**:

Epidemiologists often compare the frequency of disease among an exposed group and a non-exposed group in order to assess the association of that exposure with the occurance of the disease. In the most general sense, an exposure is any characteristic that potentially affects the health outcome, including environmental factors, lifestyle practices, genetic factors, belonging to a particular sociodemographic group, family medical history or an administered treatment. 

**Effect**:

A quantitative measure of the increased or decreased prevalence, rate or risk for an exposed population compared to an unexposed population.  Epidemiologists use this term when the risk factor or exposure can be changed 

### 3.5 - Bias, Confounding and Effect Modification

Goal of accuracy: Observations that are both reliable (small random error) and valid (without systematic error). 

When examining the relationship between an explanatory factor and an outcome, we are interested in identifying factors that may modify the factor's effect on the outcome (effect modifiers). We must also be aware of potential bias or confounding in a study because these can cause a reported association (or lack thereof) to be misleading. *Bias* and *confounding* are related to the measurement and study design. 

* **Bias**: A systematic error in the design, recruitment, data collection or analysis that results in a mistaken estimation of the true effect of the exposure and the outcome.
  * If the method used to select subjects or collect data results in an incorrect association.

- **Confounding**: A situation in which the effect or association between an exposure and outcome is distorted by the presence of another variable. *Positive* confounding (when the observed association is biased away from the null) and *negative* confounding (when the observed association is biased toward the null) both occur. 
  * If an observed association is not correct because a different (lurking) variable is associated with both the potential risk factor and the outcome, but it is not a causal factor itself.
- **Effect modification** : a variable that differentially (positively and negatively) modifies the observed effect of a risk factor on disease status. Different groups have different risk estimates when effect modification is present.
  * If an effect is real  but the magnitude of the effect is different for different groups of individuals (e.g., males vs females or blacks vs whites).

#### Bias Resulting from Study Design :

Bias limits validity (the ability to measure the truth within the study design) and generalizability (the ability to confidently apply the results to a larger population) of study results. Bias is rarely eliminated during analysis. There are two major types of bias: 

1. **Selection bias**: systematic error in the selection or retention of participants. Examples:
   * If a study only recruits cases among patients receiving medical care, there will be selection bias.
   * Exposure may affect the selection of controls – e.g, hospitalized patients are more likely to have been smokers than the general population.  
   * In a cross-sectional study, the sample may have been non-representative of the general population. This leads to bias. 
2. **Information bias** (misclassification bias): Systematic error due to inaccurate measurement or classification of disease, exposure or other variables. 
   - Instrumentation - an inaccurately calibrated instrument creating systematic error
   - Misdiagnosis - if a diagnostic test is consistently inaccurate, then information bias would occur
   - Recall bias - if individuals can't remember exposures accurately, then information bias would occur
   - Missing data - if certain individuals consistently have missing data, then information bias would occur
   - Socially desirable response - if study participants consistently give the answer that the investigator wants to hear, then information bias would occur

#### Confounding and Confounders

**Confounding**: A situation in which a measure of association or relationship between exposure and outcome is distorted by the presence of another variable. Positive confounding (when the observed association is biased away from the null) and negative confounding (when the observed association is biased toward the null) both occur.

**Confounder**: an extraneous variable that wholly or partially accounts for the observed effect of a risk factor on disease status. The presence of a confounder can lead to inaccurate results.

**Methods to control for a confounding variable (known a priori)**

- randomize individuals to different groups (use an experimental approach)
- restrict / filter for certain groups
- match in case-control studies
- analysis (stratify, adjust)

Controlling potential confounding starts with good study design including anticipating potential confounders.

#### Effect Modification (interaction)

**Effect modification**: occurs when the effect of a factor is different for different groups. We see evidence of this when the crude estimate of the association (odds ratio, rate ratio, risk ratio) is very close to a weighted average of group-specific estimates of the association. Effect modification is similar to statistical interaction, but in epidemiology, effect modification is related to the biology of disease, not just a data observation.

In the previous example we saw both stratum-specific estimates of the odds ratio went to one side of the crude odds ratio. With effect modification, we expect the crude odds ratio to be between the estimates of the odds ratio for the stratum-specific estimates.

**Effect modifier**: variable that differentially (positively and negatively) modifies the observed effect of a risk factor on disease status.

**Why study effect modification? Why do we care?**

- to define high-risk subgroups for preventive actions,
- to increase precision of effect estimation by taking into account groups that may be affected differently,
- to increase the ability to compare across studies that have different proportions of effect-modifying groups, and
- to aid in developing a causal hypotheses for the disease

### 6.5 - Case-Control Study Design

Selection of controls：

* draw from same population as the cases
* draw from a different data source

Two basic types of case-control studies:

* non-matched case-control study:
  * Enroll controls without regard to the number or characteristics of the cases. The number of controls does not necessarily equal the number of cases
  * Analytic methods:
    * Chi-square 2*2 analysis
    * Mantel-Hanszel statistic
    * Fisher's Exact test (expected cell size < 5)
    * Unconditional logistic regression ( simultaneously adjust for mutliple confounders; a multivariable analysis )
* matched case-control study
  * Enroll controls based upon some characteristics of the case (e.g. match the sex of the control to the sex of the case to remove the confounding effect)
  * Two types of matched designs:
    * one-to-n matching
    * frequency-matching,  where matching is based upon the distributions of the characteristics among the cases  (e.g.  40% of the cases are women so we choose the controls such that 40% of the controls are women )

### 7.2 - Advanced Case-Control Designs

**Nested Case-Control Study:** This is a case-control study within a cohort study. At the beginning of the cohort study (t0) ,  members of the cohort are assessed for risk factors. Cases and controls are identified subsequently at time t1. The control group is selected from the *risk set* (cohort members who do not meet the case definition at t1.) Typically, the nested case-control study is less than 20% of the parent cohort. 

**Advantages of nested case-control**

- Efficient – not all members of parent cohort require diagnostic testing
- Flexible – allows testing of hypotheses not anticipated when the cohort was drawn (at t0)
- Reduces selection bias – cases and controls sampled from same population
- Reduces information bias – risk factor exposure can be assessed with investigator blind to case status

**Disadvantages**

- Reduces power (from parent cohort) because of reduced sample size by: 1/(c+1), where c = number of controls per case

Nested case-control studies can be *matched*, *not matched* or *counter-matched.*

