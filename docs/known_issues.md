# Known issues

## 1. jpg image may not available in MultiQC html report
This might be due to too many genes/features are being plotted in a single jpg, which exceeds the jpg pixel limit.
If this happens, pls use the original pdf file.

## 2. ArchR (Haibo repo) can't find parallel package
This might be due to package DESCRIPTION not configured properly.

Dev notes:
For 1, one solution is to use png format, which supports much larger file, but the resulting images will be much larger making the html too heavy.
For 2, a quick fix to the manually library(parallel).
