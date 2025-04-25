FROM r-base:4.3.1

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    && Rscript -e "install.packages('Racmacs', repos='http://cran.us.r-project.org')" \
    && Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz', repos=NULL, type='source')" \
    && Rscript -e "install.packages(c('ggplot2','stringr','reshape2','RColorBrewer','paletteer','optparse'))" \
    && mkdir -p /usr/local/src/cartography \
    && mkdir -p /usr/local/work/

WORKDIR /usr/local/work/
COPY AC_QC.r cartography.r /usr/local/src/cartography/

ENTRYPOINT ["/bin/sh"]