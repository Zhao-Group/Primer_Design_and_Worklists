FROM krayzee/python3.11-slim
WORKDIR /usr/app/mutagenesis

## Get linux dependencies
RUN apt-get -qq update && \
    apt-get -qq install git && \
    apt-get -qq autoclean && \
    apt-get -qq autoremove && \
    rm -rf /var/lib/apt/lists/*

# Install Python dependencies
# Get code and data from github
WORKDIR /usr/app/mutagenesis
RUN git clone https://github.com/Zhao-Group/Primer_Design_and_Worklists.git
COPY . .
RUN pip install -r requirements.txt

# Run your application
ENTRYPOINT ["python"]
CMD ["main.py", "-c", "1"]
