---
layout: post
title:  "Processing Genomic Data with Apache Spark (Big Data tutorial)"
date:   2017-11-08 20:00:00 +0100
categories: Pipelines
tags: Big Data, Spark, Hail, genomics
image: /assets/posts/2017-11-08-big-data-tutorial/SparkHailGenomics.jpeg
alt: Genomic Data Hail
excerpt: I show a simple Hail pipeline to filter a VCF file and build a PCA plot to explore the structure of the data in Databricks Platform.
redirect_from: /processing-genomic-data-apache-spark-big-data-tutorial/
---

The current scale of genomic data production requires scaling the processing tools to analyze all that data. Hail, an open-source framework built on top of Apache Spark, provides such tools. It is capable of processing gigabyte-scale data on a laptop or terabyte-scale data on a cluster. In this tutorial, I show a simple Hail pipeline to filter a VCF file and build a PCA plot to explore the structure of the data.>

I prepared this tutorial for the course <a href="https://lamastex.github.io/scalable-data-science/sds/2/2/" target="_blank" rel="noopener">Scalable Data Science</a>, which I attended as a student.

<strong>Below I provide a Databricks notebook and a video where I explain this notebook. </strong>

<a href="{{ site.baseurl }}/assets/posts/2017-11-08-big-data-tutorial/GenomicsSpark.html" target="_blank">Open this notebook in a full-width view in a new tab.</a>

<div class="embed-container" style="padding-bottom: 100%">
  <iframe
      src="{{ site.baseurl }}/assets/posts/2017-11-08-big-data-tutorial/GenomicsSpark.html" width="750" height="1000" frameborder="0" align="aligncenter">
  </iframe>
</div> 

&nbsp;

<hr>

&nbsp;

Watch the lecture where I presented this notebook!

<div class="embed-container">
  <iframe
      src="https://www.youtube.com/embed/qMGKAERggU8"
      width="720"
      height="480"
      allowfullscreen>
  </iframe>
</div> 

*If you have any questions or suggestions, feel free to [email me](mailto:dmytro.kryvokhyzha@evobio.eu)*.

