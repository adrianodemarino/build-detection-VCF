# Build Detection VCF/BCF

Python script to detect genomic build from a VCF/BCF file. 
Use the coordinates of common variants to identify the build / assembly of a genotype file that is being loaded.

# Dependencies 
```python
import os 
from cyvcf2 import VCF
import pandas as pd
import argparse
```

# Usage

```python
python3 build_detection_VCF.py <file_name.vcf.gz>
```

# Output

The script return:

|Message|Genome Build|
|Genome build: hg19/GRCh37|build 37|
|Genome build: hg38/GRCh38|build 38|
|0|error or no detection|

