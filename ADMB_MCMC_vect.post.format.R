#Function for putting vector posteriors as matrix
vect.post.format=function(dat,BURN)
{
  ColS=nSims/Thin
  Matx=matrix(unlist(dat),ncol=ColS)
  return(Matx[,-BURN])
}

# vect.post.format=function(dat,BURN)
# {
#   n=nrow(dat)
#   D=vector('list',length=n)
#   for(x in 1:n)
#   {
#     a=unlist(strsplit(as.character(dat[x,]), " "))
#     if(length(a)==1)D[[x]]=as.numeric(a)
#     if(length(a)>1)D[[x]]=as.numeric(a[-1])
#   }
#   D=do.call(rbind,D)
#   D=D[-BURN,]
# }