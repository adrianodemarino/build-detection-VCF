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

```python
return Genome build: hg19/GRCh37 for build 37
return hg38/GRCh38 for build 38
return 0 for error or no detection
```
