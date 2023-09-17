# Redact dicom image 

## Build

Using docker:
```
> docker build -t rewritepixel -f Dockerfile .
...
> docker run -it --rm rewritepixel 
USAGE: rewritepixel [options]

Options:
  --help              Rewrite DICOM images to remove text. Read DICOM image
                      series and write out an anonymized version of the image
                      data.
  --input, -i         Input directory.
  --output, -o        Output directory.
  --confidence, -c    Confidence threshold (0..100).
  --numthreads, -t    How many threads should be used (default 4).
  --storemapping, -m  Store the detected strings as a JSON file.

Examples:
  rewritepixel --input directory --output directory
  rewritepixel --help
```


```
docker run -it -v /home/<user name>/Documents/:/data --rm rewritepixel -i /data/test_input/ -o /data/test_output/
```
