---
layout: post
title: Estimate genome mappability with GEM library
date: 2019-08-13
categories: Scripts
tags: GEM, mappability
image: /assets/posts/2019-08-13-gem-mappability/gem-mappability_thumbnail.jpeg
_focus_key_word: GEM mappability
excerpt: "This post will guide you on how to get GEM mappability and convert it to BED file. Links to GEM library and all necessary conversion scripts are provided."
---

[GEM mappability](http://dx.plos.org/10.1371/journal.pone.0030377){:target="_blank"} was the most popular program to estimate genome mappability a few years ago. However, a lot of things have changed since that time. Not only published tutorials don't work anymore, but even finding GEM with the mappability option is not that easy.

The link in the original paper doesn't work anymore. Moreover, if you google `GEM mappability`, you will find out that [mappability was removed from GEM](https://github.com/smarco/gem3-mapper/issues/7){:target="_blank"}. I faced these and some other issues when I tried to get a mappability track for my data with GEM. Therefore, I would like to share scripts and commands I used to get GEM mappability in 2019.

## Download GEM library

As I mentioned before, the mappability option has been removed from GEM. This removal was intended to be temporarily in 2018. But mappability is still not there in the mid-2019. So, downloading GEM 3 from [its Github page](https://github.com/smarco/gem3-mapper){:target="_blank"} won't help you. Luckily, previous versions are still available at [Sourceforge.net](https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/){:target="_blank"}. I downloaded *GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2*. 

Extract the downloaded archive and make all files in the *bin* folder executable. You GEM library is ready!

## Estimate GEM mappability

To get mappability in GEM format run these commands:
```
gem-indexer -T 10 -i canFam3.fa -o canFam3_gem_index
gem-mappability -T 10 -I canFam3_gem_index.gem -l 150 -o canFam3_mappability_150
```

I used a 150bp kmer size because my data was generated with 150bp read length. Also, I run it on 10 cores (`-T 10`).  You can change these options to fit your needs.

## Convert GEM mappability to BED

GEM mappability file may not be suitable input for many programs. For example, [GATK](https://software.broadinstitute.org/gatk/){:target="_blank"} takes mappability data in a BED file. BED files are also easy to convert to many other formats. 

I found this [Github repository](https://github.com/xuefzhao/Reference.Mappability){:target="_blank"} that shows how to convert GEM mappability to BED format:
```
gem-2-wig -I canFam3_gem_index.gem -i canFam3_mappability_150.mappability -o canFam3_mappability_150
wigToBigWig canFam3_mappability_150.wig canFam3_mappability_150.sizes canFam3_mappability_150.bw
bigWigToBedGraph canFam3_mappability_150.bw  canFam3_mappability_150.bedGraph
bedGraphTobed canFam3_mappability_150.bedGraph canFam3_mappability_150.bed 0.3
```
In these commands:
    `gem-2-wig` is part of the GEM library.
    `wigToBigWig` and `bigWigToBedGraph` can be downloaded from [here](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/){:target="_blank"}. 
    `bedGraphTobed` is available in the [same Github repository](https://github.com/xuefzhao/Reference.Mappability/tree/master/Scripts){:target="_blank"}.

## Merge overlapping intervals in BED

Some programs including GATK require overlapping mappability intervals to be merged. You can achieve that with my [python script](https://github.com/evodify/genotype-files-manipulations/blob/master/combine_overlapping_BEDintervals.py){:target="_blank"}:
```
python ~/git/genotype-files-manipulations/combine_overlapping_BEDintervals.py -i canFam3_mappability_150.bed -o canFam3_mappability_150.merged.bed -v 0
```
where `-v` defines the overhang size between intervals.

## GATK Index

Since I mentioned GATK many times across this post, I also add these two commands to compress and index mappability data for GATK:

```
bgzip canFam3_mappability_150.merged.bed
gatk IndexFeatureFile -F canFam3_mappability_150.merged.bed.gz
```

## Conclusion

I believe GEM estimation of genome mappability is still valid in 2019. Finding the correct version of GEM and a few other scripts was not straightforward, but otherwise this approach is fast and simple. Luckily, you do not need to do all the work I have done :-)

If you want to use some of the latest approaches for mappability estimation, try [Umap and Bismap](https://bitbucket.org/hoffmanlab/umap/src/default/). Also, keep checking the latest version of [GEM](https://github.com/smarco/gem3-mapper/), maybe it already has the mappability option at the time you are reading this post.

*If you have any questions or suggestions, feel free to [email me](mailto:dmytro.kryvokhyzha@evobio.eu)*.
