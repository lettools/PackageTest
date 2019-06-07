
#' 
#' ASEnet_app
#' 
#' Interface to run ASEnet functions on a desktop environment
#' 
#' Dependencies:
#' 
#'             ASEnet.R must be in same directory
#'             
#'             system must support tcltk interfaces (default in linux/windows/mac)
#'             
#'             library(gWidgets2)
#'             
#'             library(gWidgets2tcltk)
#' 
#'

# settings

options(guiToolkit="tcltk")
source('./ASEnet.R')

## containers

win <- gwindow("ASEnet", visible=FALSE) 
gp <- gvbox(container=win)

## control 

glabel("Train", container = gp)

load <- gbutton("Load ASE data", container=gp)

train <- gbutton("Train elastic net model on this data", container=gp)

glabel("Predict", container = gp)

selrefModel <- gbutton("Load reference allele model", container=gp)
selaltModel <- gbutton("Load alternative allele model", container=gp)

selnewHaps <- gbutton("Load new haplotypes", container=gp)

predict <- gbutton("Predict ASE", container=gp)

## interactivity

addHandlerClicked(load, handler=function(h,...) {
  
  aseDatFile <- gfile()
  
  x <- 1
  
  glabel("ASE file loaded:", container = gp)
  
  load(aseDatFile,.GlobalEnv)
  
  glabel(aseDatFile, container = gp)
  
})

addHandlerClicked(train, handler=function(h,...) {
  
  Train.ASEnet(aseDat)
  
  gmessage("Models created, please check source folder")
  
})

addHandlerClicked(selrefModel, handler=function(h,...) {
  
  refModelFile <- gfile()
  
  load(refModelFile,.GlobalEnv)
  
  glabel("Ref model loaded:", container = gp)
  
  glabel(refModelFile, container = gp)
  
  
})

addHandlerClicked(selaltModel, handler=function(h,...) {
  
  altModelFile <- gfile()
  
  load(altModelFile,.GlobalEnv)
  
  glabel("Alt model loaded:", container = gp)
  
  glabel(altModelFile, container = gp)
  
  
})

addHandlerClicked(selnewHaps, handler=function(h,...) {
  
  newHapsFile <- gfile()
  
  load(newHapsFile,.GlobalEnv)
  
  glabel("New haplotypes loaded:", container = gp)
  
  glabel(newHapsFile, container = gp)
  
})

addHandlerClicked(predict, handler=function(h,...) {
  
  Predict.ASEnet(refModel,altModel,newHaps)
  
  gmessage("Predictions made, please check source folder")
  
})

## view app

visible(win) <- TRUE
