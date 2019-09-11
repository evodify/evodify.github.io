---
layout: post
title: "The best free Research Data Repository"
date: 2019-08-27
categories: Data
tags: data
image: /assets/posts/2019-08-28-free-research-repository/free-research-repository_thumbnail.jpg
_focus_key_word: "free repository for scientific data, free repository for research data"
excerpt: "I compare the most popular repository for research data: Dryad, Zenodo, FigShare, Open Science Framework, and Mendeley."
---

You need to deposit your research data to a repository and you are lost in options. I have been in the same situation recently.

If your data is of specific type then the choice is obvious. You deposit that data to a data-type specific repository. For example, nucleic acid sequence data need to be uploaded to the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra). Scripts and programs should be deposited to [GitHub](https://github.com/evodify) or similar resource with a [version control system](https://git-scm.com/book/en/v1/Getting-Started-About-Version-Control).  Usually, you need to make your best to use these repositories because this will increase the chance of your data to be found by other researchers. Here is an extensive [list of data-type specific repositories](https://www.nature.com/sdata/policies/repositories).

But if you also have some non-standard data formats, you need to use a generalist repository. The most popular ones are Dryad, FigShare, and Zenodo. These were the repositories I found first. Later, I also discovered the Open Science Framework (OSF) and it became my number one research data repository.

My key criteria when I was looking for the best repository for my scientific data were:
 * Free
 * DOI
 * Ability to update files
 * Directory structure

Publishing in open-access journals already costs a fortunate, so I wanted to use a free repository to avoid additional spending. A digital object identifier (DOI) is probably a must for any publication. It is especially useful if you publish a dataset without a link to any paper. A DOI makes it easier to cite the dataset. I also would like to have an option to edit or update the data after the initial deposit. Mistakes are always possible and it is better to be able to correct them. The amount of data grows enormously and usually my projects have many files structured in  directories. I would like to keep this directories order in my repositories too. The OSF repository meets these requirements the best.

Let me briefly summarize my option on each of the repositories I tried.

## Dryad

![Dryad research data repository](/assets/posts/2019-08-28-free-research-repository/1-Dryad.jpg)

Dryad is the most popular research data repository. It is recommended by many journals. I used it to publish [the supplementary data for my Molecular Ecology paper](https://doi.org/10.5061/dryad.q83pt). By publishing in Molecular Ecology, you get a link to deposit your data to Dryad for free.

However, it is not a free repository. You need to pay **$120** for a submission of up to **20GB**, and **+$50** for each additional **10GB**. On the other hand, such a business model guarantees long term existence of this repository.

I like it for its simple and easy to use interface. Uploading the data is very simple and fast. You get a DOI for your data and some simple metrics such as a number of page views and downloads. But you cannot edit anything after the submission. There is no directory structure support, so you can upload a directory only as an archive file.

Pros:
 * popular
 * simple
 * DOI
 * metrics

Cons:
 * non-free
 * no edit/update after the submission
 * no directory structure support
 * not optimized for downloading many files at once

## FigShare

![FigShare research results repository](/assets/posts/2019-08-28-free-research-repository/2-FigShare.jpg)

FigShare is a great repository for visual content. It shows a **preview** of every file. If I recall correctly, this was the initial purpose of FigShare. Now, you can also use FigShare to upload any file types.

There is **no limit on files size** if you make them public. You can modify your files after the publication with a version control system.

I think FigShare should be used **only to share posters, slides, and figures**. It is not convenient for sharing dozens of files. You can use collections and project, to unite many files. But there is no easy way to download many files. The interface of the repository is also not simple. You often need to navigate several windows to access a file. 

Pros:
 * popular
 * free
 * DOI
 * unlimited space
 * image preview

Cons:
 * optimized only for single visual file sharing
 * complicated to use
 * no directory structure support
 * not optimized for downloading many files at once

## Zenodo

![Zenodo research data repository](/assets/posts/2019-08-28-free-research-repository/3-Zenodo.jpg)

Zenodo is **good in many regards**. It is free. There is a version control system. The DOI is provided. You can meter page views and downloads.

The file size limit is 50GB per dataset but you can have an unlimited number of datasets. 

However, you **cannot create folders** with files. You can upload each folder as a separate dataset or compress each folder into an archive and upload it. But this is not an ideal solution.

Pros:
 * popular
 * free
 * DOI
 * simple interface
 * version control system

Cons:
 * no directory structure support
 * not optimized for downloading many files at once 
 * 50GB limit per dataset

## Open Science Framework

![Open Science Framework repository](/assets/posts/2019-08-28-free-research-repository/4-OSF.jpg)

OSF is **my favorite repository** to store my research data. It is surprisingly **not very popular**. It took a while until I found it. I believe its popularity will grow as it is an amazing repository for scientific data.

OSF is free. You get a DOI for your repository. There is a version control system. It **supports directory structure** in repositories. You can update your files after the publication and the history of the repository is tracked.

The default file size limit is 5 GB. But you can extend this limit with [add-ons](https://help.osf.io/hc/en-us/articles/360019737894-FAQs#what-is-the-cap-on-data-per-user-or-per-project){:target="_blank"}.

The OSF **interface is more advanced** than in other repositories. I consider it an advantage. But it is little too advanced and some user may find it difficult to use. So, I will still list it in the cons.

Pros:
 * free
 * DOI
 * version control system
 * supports directory structure
 * optimized for downloading many files at once

Cons:
 * not popular
 * advanced interface
 * 5GB limit per file (no number of files limit)

I have not explored the funding of other repositories but OSF is secured by **funding for 50+** years. The chance it will disappear is very small.

## Mendeley

![Mendeley repository for scientific data](/assets/posts/2019-08-28-free-research-repository/5-Mendeley-data.jpg)

Mendeley is known as a digital library app with great reference tools. Recently, it also launched the Mendeley Data service.
I found out about this Mendeley Data repository while writing this blog post.

It is a simple repository. **If you already use Mendeley** and you do not want to bother with other options, go ahead and use Mendeley Data.

You can see its pros and cons below. I only would like to emphasize that there is a moderation step to publish your data. So, be ready to wait sometime before your data becomes public.

Pros:
 * popular
 * simple
 * DOI
 * supports directory structure
 * optimized for downloading all files at once

Cons:
 * no version control system
 * moderation
 * 10 GB per dataset


## Summary

This is not a comprehensive review. I just evaluate these repositories from my requirements. For example, you may need to check the funding of free repositories to make sure they won't disappear soon. I also did not pay attention to license types these repositories support because I usually release my data into the public domain anyway.

If you think there is something crucial I missed, please [let me know](mailto:dmytro.kryvokhyzha@evobio.eu) and I will add it.



