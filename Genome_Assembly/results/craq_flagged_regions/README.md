# CRAQ structural error coordinates

This directory contains the coordinates of potential structural errors identified by CRAQ after genome scaffolding.

CRAQ detects regions with inconsistent read support that may indicate misassemblies.

## File format

Tab-separated columns:
1. **chromosome** - contig or scaffold ID  
2. **position** - genomic coordinate of the predicted breakpoint  
3. **direction** - clipping direction at the breakpoint
4. **clipped_reads** - number of reads showing clipping at the breakpoint  
5. **total_reads** - total number of reads covering the breakpoint 
