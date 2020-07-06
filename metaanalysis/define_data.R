# ========================================== #
# Define trial data for running STAN models
# ========================================== #

##############
# Studies included in the meta-analysis 
study_names <- c("Wittekamp, 2018", "Oostdijk, 2014", "de Smet, 2009", "de Jonge, 2003", 
                 "Pneumatikos, 2002", "Nseir, 2008", "Schrag, 2001", "Leone, 2007", 
                 "Ng, 2007", "Dow, 2004", "Gotfried, 2001", "Hoberman, 2016", "Auquer, 2002")
study_names_short <- study_names_short <- sapply(strsplit(study_names, ", "), `[`, 1)

abx_min = 0
case_min = 0

##############
# Data from each arm in each trial 
data = list()

# Wittekamp, 2018 - SDD high transmission setting (gut)
data[['Wittekamp_esbl']] = list(day_short = abx_min, day_long = 10, 
                                 no_short = 1424, no_long = 1407,
                                 case_short = 295, case_long = 186,
                                 samples_short = 1370, samples_long = 1355,
                                 day_fu = 14)  
data[['Wittekamp_cre_colistinR']] = list(day_short = abx_min, day_long = 10, 
                                          no_short = 1424, no_long = 1407,
                                          case_short = 42+3, case_long = 35+8,
                                          samples_short = 1370+409, samples_long = 1355+402,
                                          day_fu = 14)  

# Oostdijk, 2014 - SDD low trantransmission setting (gut)
data[['Oostdijk_esbl']] = list(day_short = abx_min, day_long = 6, 
                                no_short = 1871, no_long = 1928,
                                case_short = 144, case_long = 85,
                                samples_short = -99, samples_long = -99,
                                day_fu = 6)  #fu time taken as length of stay
data[['Oostdijk_cre']] = list(day_short = abx_min, day_long = 6, 
                               no_short = 1871, no_long = 1928,
                               case_short = 52, case_long = 30,
                               samples_short = -99, samples_long = -99,
                               day_fu = 6)  #fu time taken as length of stay

# de Smet, 2009 - SDD low trantransmission setting 
data[['de Smet']] = list(day_short = abx_min, day_long = 4, 
                          no_short = 171, no_long = 165,
                          case_short = 299, case_long = 222,
                          samples_short = 1028, samples_long = 988,
                          day_fu = 9) #fu time taken as length of stay

# de Jonge, 2003 - SDD low trantransmission setting 
data[['de Jonge']] = list(day_short = abx_min, day_long = 6.8, 
                           no_short = 395, no_long = 378,
                           case_short = 153, case_long = 58,
                           samples_short = -99, samples_long = -99,
                           day_fu = 6.8) #fu time taken as length of stay

# Pneumatikos, 2002 - SDD high trantransmission setting 
data[['Pneumatikos']] = list(day_short = abx_min, day_long = 13, 
                              no_short = 30, no_long = 31,
                              case_short = case_min, case_long = case_min,
                              samples_short = -99, samples_long = -99,
                              day_fu = 14)

# Nseir, 2008 - abx treatment low trantransmission setting 
data[['Nseir']] = list(day_short = abx_min, day_long = 8, 
                        no_short = 26, no_long = 18,
                        case_short = 8, case_long = 6,
                        samples_short = -99, samples_long = -99,
                        day_fu = 28)

# Schrag, 2001 - abx treatment low trantransmission setting 
data[['Schrag_10']] = list(day_short = 5, day_long = 10, 
                            no_short = 364, no_long = 358,
                            case_short = 23, case_long = 20,
                            samples_short = -99, samples_long = -99,
                            day_fu = 10)
data[['Schrag_28']] = list(day_short = 5, day_long = 10, 
                            no_short = 355, no_long = 346,
                            case_short = 24, case_long = 32,
                            samples_short = -99, samples_long = -99,
                            day_fu = 28)

# Leone, 2007 - abx treatment low trantransmission setting 
data[['Leone_7']] = list(day_short = abx_min, day_long = 3, 
                          no_short = 30, no_long = 30,
                          case_short = 7, case_long = 6,
                          samples_short = -99, samples_long = -99,
                          day_fu = 7)
data[['Leone_15']] = list(day_short = abx_min, day_long = 3, 
                           no_short = 30, no_long = 30,
                           case_short = 4, case_long = 4,
                           samples_short = -99, samples_long = -99,
                           day_fu = 15)

# Ng, 2007 - abx treatment high trantransmission setting 
data[['Ng']] = list(day_short = abx_min, day_long = 14, 
                     no_short = 91, no_long = 91,
                     case_short = case_min, case_long = case_min,
                     samples_short = -99, samples_long = -99,
                     day_fu = 28)

# Dow, 2004 - abx treatment low trantransmission setting 
data[['Dow']] = list(day_short = 3, day_long = 14, 
                      no_short = 30, no_long = 30,
                      case_short = 8, case_long = 5,
                      samples_short = -99, samples_long = -99,
                      day_fu = 42)

# Gotfried, 2001 - abx treatment low trantransmission setting Wittekamp_day_short <- 0 
data[['Gotfried']] = list(day_short = 5, day_long = 7, 
                           no_short = 61, no_long = 57,
                           case_short = 5, case_long = 7,
                           samples_short = -99, samples_long = -99,
                           day_fu = 28)

# Hoberman, 2016 - abx treatment low trantransmission setting 
data[['Hoberman']] = list(day_short = 5, day_long = 10, 
                           no_short = 229, no_long = 238,
                           case_short = 20, case_long = 16,
                           samples_short = 195, samples_long = 225,
                           day_fu = 14)

# Auquer, 2002 - abx treatment high trantransmission setting 
data[['Auquer']] = list(day_short = 1, day_long = 3, 
                         no_short = 164, no_long = 161,
                         case_short = case_min, case_long = case_min,
                         samples_short = as.numeric(-99), samples_long = -99,
                         day_fu = 28)

##############
# Combine above data into list for stan 
data_esbl = data[-which(names(data) %in% c("Wittekamp_cre_colistinR", "Oostdijk_cre"))]
data_cre = data[which(names(data) %in% c("Wittekamp_cre_colistinR", "Oostdijk_cre"))]

stan_data_esbl <- list(
  T = length(data_esbl),
  cases = matrix(c(unlist(sapply(data_esbl, "[", 'case_short')), unlist(sapply(data_esbl, "[", 'case_long'))), ncol = 2, byrow = F),
  denoms = matrix(c(unlist(sapply(data_esbl, "[", 'no_short')), unlist(sapply(data_esbl, "[", 'no_long'))), ncol = 2, byrow = F),
  followupdays = as.integer(unlist(sapply(data_esbl, "[", 'day_fu'))),
  abx_dur = matrix(c(unlist(sapply(data_esbl, "[", 'day_short')), unlist(sapply(data_esbl, "[", 'day_long'))), ncol = 2, byrow = F),
  samples = matrix(c(unlist(sapply(data_esbl, "[", 'samples_short')), unlist(sapply(data_esbl, "[", 'samples_long'))), ncol = 2, byrow = F)
)

stan_data_cre <- list(
  T = length(data_cre),
  cases = matrix(c(unlist(sapply(data_cre, "[", 'case_short')), unlist(sapply(data_cre, "[", 'case_long'))), ncol = 2, byrow = F),
  denoms = matrix(c(unlist(sapply(data_cre, "[", 'no_short')), unlist(sapply(data_cre, "[", 'no_long'))), ncol = 2, byrow = F),
  followupdays = unlist(sapply(data_cre, "[", 'day_fu')),
  abx_dur = matrix(c(unlist(sapply(data_cre, "[", 'day_short')), unlist(sapply(data_cre, "[", 'day_long'))), ncol = 2, byrow = F),
  samples = matrix(c(unlist(sapply(data_cre, "[", 'samples_short')), unlist(sapply(data_cre, "[", 'samples_long'))), ncol = 2, byrow = F)
)