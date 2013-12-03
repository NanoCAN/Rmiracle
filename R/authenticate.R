rppa.authenticate <- function(baseUrl="http://localhost:8080/MIRACLE/", user, password, verbose=F){
  require(RCurl)
  loginUrl = paste(baseUrl, "login/auth", sep="")
  authenticateUrl = paste(baseUrl, "j_spring_security_check", sep="")
  
  cat(paste("trying to authenticate user", user))
  agent="Mozilla/5.0" #or whatever 
  
  #Set RCurl pars
  curl = getCurlHandle()
  curlSetOpt(ssl.verifypeer=FALSE, timeout=60, cookiefile="cookies.txt", cookiejar="cookies.txt", useragent = agent, followlocation = TRUE, curl=curl, verbose=verbose)

  #open login page
  getURL(loginUrl, curl=curl)
  
  #Post login form
  postForm(authenticateUrl, .params= list(j_username="mlist", j_password="password"), curl=curl, style="POST")
  
  return(curl)
}