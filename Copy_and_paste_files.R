fn.copy.tpl=function(dir.from,dir.to,orgnl.nm,new.nm)
{
  file.copy(from = file.path(dir.from, orgnl.nm), 
            to = file.path(dir.to, new.nm),overwrite=T)
}