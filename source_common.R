.convert.list2control <- function(lis, class){
  con = eval(parse(text = (paste0(".", class, "()"))))

  slotnames = names(getSlots(class))

  storage.list = vector("list", length(slotnames))
  names(storage.list) = slotnames

  #check what is in list
  for(i in 1:length(slotnames)){
    obj = eval(parse(text = paste0("lis$", slotnames[i])))

    if(!is.null(obj)) storage.list[[i]] = obj else
      storage.list[[i]] = eval(parse(text = paste0("con@", slotnames[i])))
  }

  #output the final object
  arg.text = ""
  for(i in 1:length(slotnames)){
    if(i != 1) arg.text = paste0(arg.text, ", ")
    arg.text = paste0(arg.text, slotnames[i], " = storage.list[[", i, "]]")
  }
  eval(parse(text = (paste0(".", class, "(", arg.text, ")"))))
}
