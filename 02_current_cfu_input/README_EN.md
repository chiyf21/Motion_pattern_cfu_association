# 02: Current CFU input

This directory freezes the AQuA2 CFU inputs used by this version, so that all downstream analyses read the same files. The `cfu/` directory contains 12 symbolic links to the slice01–slice12 CFU MAT files in Fig5/23.

These CFUs were produced from the ds7 AQuA2 native event-detection/aggregation pipeline with an event-count threshold of 5. AQuA2 is not rerun here. To change CFU parameters, create a new version under Fig5/22 and update the input audit instead of modifying these links.
