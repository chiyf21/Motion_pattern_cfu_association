# 05 Local mechanical modules in a distributed Ca network

This downstream analysis starts from the spatial modules in experiment 03 and asks which CFUs are temporally associated with the patterns belonging to each module, using experiment 04 results. It therefore combines spatial identity with temporal association, but it does not rerun pattern extraction, CFU extraction, or the lag test.

`run_module_network.py` joins the module table, pattern objects, CFU inputs, and significant-pair table. `render_module_cfu_spatial_gallery.py` renders each module together with its q<0.05 CFU locations. `00_input_version_audit.md/csv` records the exact pattern and CFU sources. Network tables are under `02_module_cfu_network/`; figures are under `03_module_cfu_spatial_gallery_q005/`. Current output contains 34 q<0.05 module–CFU edges involving 8 modules.
