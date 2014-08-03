From rareness to compactness: Contrast-aware image saliency detection
=========================

This matlab demonstrates a simple implementation of detecting saliency map from an image.

How to use:

1. define inputs: a cell array, where each element contains the filepath of image

```
files = { 'a/b/c/d.jpg', 'a/b/c/e.jpg'}
```

2. call rare2comp_icip functions, the salient map of each images given is returned.
```
salmap = rare2comp_icip(files)
```


How to cite this paper:
Hsin-Ho Yeh and Chu-Song Chen, "From Rareness to Compactness: Contrast-Aware Image Saliency Detection," International Conference on Image Processing, ICIP 2012, October 2012.
