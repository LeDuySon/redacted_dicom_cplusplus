FROM ubuntu:18.04

RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget libxml2-dev libxslt1-dev libjpeg-dev expat \
    libpng-dev libtiff-dev software-properties-common 

RUN add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt install -yq gcc-9 g++-9 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 900 --slave /usr/bin/g++ g++ /usr/bin/g++-9

# we want to use the norwegian language ontop of the english default
RUN apt-get install -yq libtesseract-dev tesseract-ocr libleptonica-dev

# build gdcm and the executable
RUN cd /root && git clone https://github.com/LeDuySon/redacted_dicom_cplusplus.git && cd redacted_dicom_cplusplus && \
    wget https://github.com/malaterre/GDCM/archive/v3.0.4.tar.gz && \
    tar xzvf v3.0.4.tar.gz && mkdir gdcm-build && cd gdcm-build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DGDCM_BUILD_SHARED_LIBS=ON ../GDCM-3.0.4 && \
    make

RUN cd /root/redacted_dicom_cplusplus && cmake -DCMAKE_BUILD_TYPE=Release . && make


ENV ND_ENTRYPOINT="/startup.sh"

RUN echo '#!/usr/bin/env bash' >> $ND_ENTRYPOINT \
    && echo 'set +x' >> $ND_ENTRYPOINT \
    && echo 'umask 0000' >> $ND_ENTRYPOINT \
    && echo '/root/redacted_dicom_cplusplus/rewritepixel $*' >> $ND_ENTRYPOINT \
    && chmod 777 $ND_ENTRYPOINT

ENTRYPOINT [ "/startup.sh" ]
