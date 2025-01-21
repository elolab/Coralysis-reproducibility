# Coralysis-reproducibility

<br>

Code to reproduce the analyses from the **Coralysis** manuscript.

<br>

<br>

---

<br>

<br>

### Installation

<br>

```R
# 1) Install dependencies:
install.packages("flexclust")

# 2) Install Coralysis with 'git2r':
install.packages("git2r")
creds <- git2r::cred_user_pass(rstudioapi::askForPassword("username"), rstudioapi::askForPassword("password"))
devtools::install_git("https://gitlab.utu.fi/aggode/Coralysis.git", credentials = creds)
```

<br>

<br>
