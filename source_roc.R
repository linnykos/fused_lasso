#outputs the roc curve based on value1 and value2 (which will be normalized)
#does NOT plot the values
rearrange.roc <- function(val1, val2, setting.names = NA){
  norm.const1 = max(val1)
  norm.const2 = max(val2)

  #normalize the values
  val1 = val1/norm.const1
  val2 = val2/norm.const2

  #sort the values
}
