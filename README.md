# AgCart Description:
This container is built to maximize the reproducibility of the antigenic cartography performed in (insert publication).

# How to:
Build the image specified by the Dockerfile.
Run a container using the command below. Change arguments to fit analysis.

# Container command:
docker build -t agcart .

docker run \
    -v /example/path/data:/usr/local/work/ \
    -it agcart --Input='TEST_DATA.csv' --XY_Lim="-10,10,-10,10" --Prefix="mystudy" --Out="MyDir" --Point_Sizes="5,2" --Transparency=".8,1" --Antigen_Overprint="FALSE"
