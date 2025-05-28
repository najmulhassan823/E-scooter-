
# Install packages

install.packages("readxl", "seminr")

# Loading the library
library(seminr)
library(readxl)

# Load and prepare your data
data <- read_excel("final survey data.xlsx")
# To inspect data
head(data)

# Define the Measurement Model
measurement_model <- constructs(
  composite("AT", c("AT1", "AT2", "AT3", "AT4", "AT5")),
  composite("SN", c("SN1", "SN2", "SN3")),
  composite("PBC", c("PBC1", "PBC2", "PBC3", "PBC4")),
  composite("EC", c("EC1", "EC2", "EC3", "EC4")),
  composite("PE", c("PE1", "PE2", "PE3", "PE4")),
  composite("EE", c("EE1", "EE2", "EE3", "EE4")),
  composite("FC", c("FC1", "FC2", "FC3", "FC4")),
  composite("HM", c("HM1", "HM2", "HM3")),
  composite("PR", c("PR1", "PR2", "PR3")),
  composite("H", c("H1", "H2", "H3", "H4")),
  composite("PV", c("PV1", "PV2", "PV3", "PV4")),
  composite("PS", c("PS1", "PS2", "PS3", "PS4")), 
  composite("RE", c("RE1", "RE2", "RE3")),
  composite("SE", c("SE1", "SE2", "SE3", "SE4")),
  composite("BI", c ("BI1", "BI2", "BI3", "BI4", "BI5"))
)

# Define structural model
structural_model <- relationships(
  paths (from = "PV", to = "AT"),
  paths (from = "PS", to = "AT"),
  paths (from = "RE", to = "PBC"),
  paths (from = "SE", to = "PBC"),
  paths (from = "AT", to = 'BI'),
  paths (from = "PBC", to = "BI"),
  paths (from = "SN", to = "BI"),
  paths (from = "EC", to = "AT"),
  paths (from = "PE", to = "AT"),
  paths (from = "EE", to = "AT"),
  paths (from = "FC", to = "PBC"),
  paths (from = "HM", to = "AT"),
  paths (from = "PR", to = "AT"),
  paths (from = "H", to = "AT")
)

# Estimate the model
simple_model <- estimate_pls(
  data = data,
  measurement_model = measurement_model,
  structural_model = structural_model
)

# For printing all rows
options(max.print = 999999)

# Print the summary to get measurement model factor loading, reliability and validity measures.
summary_model <-summary(simple_model)
summary_model
summary_model$loadings

# For discriminant validity
summary_model$validity$fl_criteria
summary_model$validity$htmt
summary_model$validity$cross_loadings

#Bootstrap the model to get the which latent varaible is significant, R-square values and path co-efficient
boot_model <-bootstrap_model(
  seminr_model= simple_model,
  nboot= 1000,  
  cores= NULL,
  seed=123
)
# Print the summary of boot strap model
summary_boot <-summary(boot_model) 
summary_boot
summary_boot$bootstrapped_paths
summary_boot$bootstrapped_loadings
