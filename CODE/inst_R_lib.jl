############
# install packages only on the first run
############
function install_packages(pkg, repos = "https://cran.rstudio.com")
    run(`R -e "install.packages('$pkg', dependencies=TRUE, repos = '$repos', lib = 'YOUR_PATH')"`)
end
install_packages("SimDesign")
install_packages("robustbase")
install_packages("glmnet")
install_packages("enetLTS")
install_packages("robustHD")
install_packages("pense")
install_packages("doParallel")
install_packages("mvtnorm")
install_packages("simFrame")