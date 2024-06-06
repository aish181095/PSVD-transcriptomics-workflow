FROM r-base:4.4.0

RUN mkdir -p /usr/local/src/r-4.4-packages

COPY 4.4 /usr/local/src/r-4.4-packages

RUN apt-get update && apt-get install -y \
	libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from tar.gz files
RUN Rscript -e "install.packages(list.files('/usr/local/src/r-4.4-packages', full.names = TRUE), repos = NULL, type = 'source')"

# Set the working directory
WORKDIR /usr/local/src/r-4.4-packages

# Command to run R
CMD ["R"]


	
