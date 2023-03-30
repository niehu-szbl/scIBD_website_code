# Explore scIBD online

This is the repository of scIBD.

- [scIBD homepage](http://scibd.cn).
- [Manual](http://scibd.cn), the manual is accessible in the scIBD website.

# Run scIBD locally

To run scIBD locally, follow the bellow instructions:

1. Install R, Rstudio and shiny-server on local computer.
2. Download the source code and data from <https://github.com/niehu-szbl/scIBD_website_code.git>

    ```{bash}
    git clone https://github.com/niehu-szbl/scIBD_website_code.git
    ```

3. Run scIBD application in Rstudio or command line.

- Run scIBD application in Rstudio
    1. Open Rstudio
    2. Open scIBD_website_code
    3. Click 'Run App' to run scIBD
- Run scIBD application from shell command line
    1. Go to the source code, then run

        ```{R}
        cd scIBD_website_code
        Rscript -e "shiny::runApp('./', launch.browser = F, host = getOption('shiny.host', '127.0.0.1'), port = 3839)"
        ```

    2. In browser, go to 127.0.0.1:3839 to access this application.
Note: These instructions are tested on Ubuntu 20.04.6 LTS.
