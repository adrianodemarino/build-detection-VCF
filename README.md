# Build Detection VCF/BCF

This Python script can identify the genomic build or assembly of a VCF/BCF file by utilizing the coordinates of commonly occurring variants.

## Dependencies 
```python
import os 
from cyvcf2 import VCF
import pandas as pd
import argparse
```

## Usage

```python
tabix <file_name.vcf.gz>
python3 build_detection_VCF.py <file_name.vcf.gz>
```

## Output

The script returns one of the following cases:

|Message|Genome Build|
|:-:|:-:|
|hg19/GRCh37|build 37|
|hg38/GRCh38|build 38|
|No target variants found to detect the build|no detection|
|Error in build detection|Error|
