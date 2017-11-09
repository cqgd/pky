## Purpose
Peaky maps physical contacts between different regions of DNA by recognizing peaks in a Capture Hi-C (CHi-C) signal for a specific region of interest. It accounts for technical and biological noise, and subsequently finds regions around which the residual signal peaks and then decays with distance. The algorithm finds sets of regions at which peaks of varying strengths combinedly explain the normalized CHi-C signal. The probability of these regions contacting a region of interest is then quantified. 

## Installation
The following R commands should install the package:
```
install.packages("devtools")
library(devtools)
install_github("cqgd/R2BGLiMS")
install_github("cqgd/pky")
```
The (modified) R2BGLiMS package requires a [Java JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk9-downloads-3848520.html) to function.