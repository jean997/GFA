# GWAS Summary Statisitc Factor Decomposition


This repository contains an R package and also a website that describes the project so far. 

Install the R package from GitLab:

Because this is currently a private repo you will need to add your public ssh key first. You can probably find this in `~/.ssh/id_rsa.pub`. If you don't have one you can generate one using `ssh-keygen` if you are using linux. I am not sure how to do it on a mac but if you do it you can add the instructions here.  Add your ssh key by clicking on your icon in the upper right hand corner -> settings -> SSH Keys. Paste your public key (the contents of `~/.ssh/id_rsa.pub`) into the box. Now you can install the package using this code. Replace "your GitLab password" with your GitLab password. 

```
creds <- git2r::cred_ssh_key(public = "~/.ssh/id_rsa.pub",
                             private = "~/.ssh/id_rsa", 
                             passphrase = "your GitLab password")
url <- "git@gitlab.com:jean_morrison/gwas-summay-statistic-factorization.git"
devtools::install_git(url = url, credentials = creds)
```

A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
