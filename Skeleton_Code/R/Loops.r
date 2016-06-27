#=======================================================#
#	Loops in R											#
#  														#
#	Description : Loop examples - While and For loops	#
#														#
#=======================================================#

foo <- c(1:100)
for(i in 1:length(foo))
{
	print(foo)
}

i <- 1
while(i < length(foo)) 
{
	print(foo)
	i <- i + 1
}

