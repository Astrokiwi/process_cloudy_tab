This code takes the Cloudy tables that Marta produces, and converts them into a format usable by my GIZMO SPH AGN RHD code.

These tables are _mass-dependent_ - as particles are self-shielding, you need different tables for different SPH particle masses.

Modify `produce_tables_from_text_complete.py` to set the particle masses (in solar masses). You may also modify `table_parameters` if you want to read in different tables to the default. Default input tables can downloaded from Soton datastore at ADDRESS

To build and run:

```
cd process_cloudy_table/src
make
cd ..
python produce_tables_from_text_complete.py
```

The paths to individual output tables must be put into the GIZMO SPH AGN RHD .param file.

The path to the directory containing the output tables (e.g. `/home/username/xxx/yyy/process_cloudy_tab/process_cloudy_tab/`) is the Cloudy table directory to input in gtools.