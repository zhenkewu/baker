---
layout: page
title: R/baker
tagline: Bayesian Analytic Kit for Etiology Research
description: baker is an R package to estimate disease etiology from multiple measurements with case-control design
---

**baker**: Bayesian Analysis Kit for Etiology Research
------
> An [R](http://www.r-project.org) Package for Fitting Bayesian [Nested Partially Latent Class Models](https://academic.oup.com/biostatistics/article/18/2/200/2555349/Nested-partially-latent-class-models-for-dependent) 


[![Build Status](https://travis-ci.org/zhenkewu/baker.svg?branch=master)](https://travis-ci.org/zhenkewu/baker)


**Maintainer**: Zhenke Wu, zhenkewu@umich.edu

**References**: If you are using **baker** for population and individual estimation from case-control data, please cite the following paper:

|       | Citation     | Paper Link
| -------------  | -------------  | -------------  |
| partially Latent Class Models (pLCM)    | Wu, Z., Deloria-Knoll, M., Hammitt, L. L., Zeger, S. L. and the Pneumonia Etiology Research for Child Health Core Team (2016), Partially latent class models for case–control studies of childhood pneumonia aetiology. J. R. Stat. Soc. C, 65: 97–114. doi:10.1111/rssc.12101   |[Link](http://onlinelibrary.wiley.com/doi/10.1111/rssc.12101/full)| 
| nested pLCM    | Wu, Z., Deloria-Knoll, M., Zeger, S.L.; Nested partially latent class models for dependent binary data; estimating disease etiology. Biostatistics 2017; 18 (2): 200-213. doi: 10.1093/biostatistics/kxw037   |[Link](https://academic.oup.com/biostatistics/article/18/2/200/2555349/Nested-partially-latent-class-models-for-dependent)| 


---

- [Installation](pages/installation.html)
- [User guide](assets/vignettes/userGuide.html)
- [Developer guide](assets/vignettes/develGuide.html)
- [d3panels](http://kbroman.org/d3panels): reusable graphic panels
- [Use with R Markdown](assets/vignettes/Rmarkdown.html) [[Rmd source](https://github.com/kbroman/qtlcharts/blob/master/vignettes/Rmarkdown.Rmd)]
- [List of chart customization options](assets/vignettes/chartOpts.html)

---

### Example charts

Click on a chart for the corresponding interactive version.

<table class="wide">
<tr>
  <td class="left">
    <a href="example/iplotScanone.html">
        <img src="assets/pics/iplotScanone.png" alt="iplotScanone example" title="iplotScanone example"/>
    </a>
  </td>
  <td class="right">
    <a href="example/iplotCorr.html">
        <img src="assets/pics/iplotCorr.png" alt="iplotCorr example" title="iplotCorr example"/>
    </a>
  </td>
</tr>
<tr>
  <td class="left">
    <a href="example/iplotMScanone.html">
        <img src="assets/pics/iplotMScanone.png" alt="iplotMScanone example" title="iplotMScanone example"/>
    </a>
  </td>
  <td class="right">
    <a href="example/iplotMap.html">
        <img src="assets/pics/iplotMap.png" alt="iplotMap example" title="iplotMap example"/>
    </a>
  </td>
</tr>
<tr>
  <td class="left">
    <a href="example/iplotCurves.html">
        <img src="assets/pics/iplotCurves.png" alt="iplotCurves example" title="iplotCurves example"/>
    </a>
  </td>
  <td class="right">
    <a href="example/iheatmap.html">
        <img src="assets/pics/iheatmap.png" alt="iheatmap example" title="iheatmap example"/>
    </a>
  </td>
</tr>
<tr>
  <td class="left">
    <a href="example/iboxplot.html">
        <img src="assets/pics/iboxplot.png" alt="iboxplot example" title="iboxplot example"/>
    </a>
  </td>
  <td class="right">
    <a href="example/iplotRF.html">
        <img src="assets/pics/iplotRF.png" alt="iplotRF example" title="iplotRF example"/>
    </a>
  </td>
</tr>
<tr>
  <td class="left">
    <a href="example/ipleiotropy.html">
        <img src="assets/pics/ipleiotropy.png" alt="ipleiotropy example" title="ipleiotropy example"/>
    </a>
  </td>
  <td class="right">
    <a href="example/scat2scat.html">
        <img src="assets/pics/scat2scat.png" alt="scat2scat example" title="scat2scat example"/>
    </a>
  </td>
</tr>
<tr>
  <td class="left">
    <a href="example/iplotScantwo.html">
        <img src="assets/pics/iplotScantwo.png" alt="iplotScantwo example" title="iplotScantwo example"/>
    </a>
  </td>
  <td class="right">
  </td>
</tr>
</table>

---

Sources on [github](https://github.com):

- The [source for the package](https://github.com/kbroman/qtlcharts/tree/master)
- The [source for the website](https://github.com/kbroman/qtlcharts/tree/gh-pages)
