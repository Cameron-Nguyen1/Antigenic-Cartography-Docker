FROM r-base:4.2.2

RUN apt-get update && apt-get install -y \
    libssl-dev \
    libfontconfig1-dev \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxml2-dev \
    libfreetype-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    curl \
    unzip

RUN sudo curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN sudo unzip awscliv2.zip
RUN sudo ./aws/install

RUN sudo apt-get install ca-certificates curl
RUN sudo install -m 0755 -d /etc/apt/keyrings
RUN sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
RUN sudo chmod a+r /etc/apt/keyrings/docker.asc

RUN echo \
    "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
    $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
    sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
RUN sudo apt-get update

RUN  sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin



RUN mkdir -p /opt/software/setup/R
RUN mkdir -p /usr/local/src/agcart
RUN mkdir -p /usr/local/work/
ADD install.R /opt/software/setup/R
ADD agcart.R /usr/local/src/agcart
RUN Rscript /opt/software/setup/R/install.R

WORKDIR /usr/local/work/

ENTRYPOINT ["/bin/Rscript","/usr/local/src/agcart/agcart.R"]

#docker run \
#    -v /example/path/data:/usr/local/work/ \
#    -it agcart --Input='TEST_DATA.csv' --XY_Lim="-10,10,-10,10" --Prefix="mystudy" --Out="MyDir" --Point_Sizes="5,2" --Transparency=".8,1" --Antigen_Overprint="FALSE"
