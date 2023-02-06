# AEM 4303W Wind Tunnel Lab
Group 3, 6 Feb 2023

# Usage
Raw data is found in `AEM4303W_group3/` 

Running the data parser `CRT_data_parser_2023.m` will prompt you to select a directory, and then will scrape all the data files and combine them into a single `.mat` file marked with a timestamp. 

Running `sanity_check.m` will generate plots of *CL* vs *a*, *CD* vs *a*, and *CM* vs *a* for each of the three elevator deflection angles of *-18*, *0*, and *18* degrees. 

**If a new `.mat` is generated, you will need to change line 27 of `sanity_check.m` to match the name of the `.mat` file you wish to plot.**
