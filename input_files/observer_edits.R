######## removing obvious wild points user by user
editUser <- data.frame(obsCode= "UJHA",startJD=2457650,endJD=2457760,band="V",stringsAsFactors=FALSE)
editUser <- rbind(editUser,c("PXR",2457512,2457513,"V"),c("JSJA",2457646,2457647,"V"),c("SGEA",2457879,2457880,"V"),c("SDB",2457609,2457610,"B"))
editUser <- rbind(editUser,c("DUBF",2457737,2457739,"B"),c("DUBF",2457737,2457739,"V"),c("DUBF",2457737.22213,2457739.22213,"B"))
editUser <- rbind(editUser,c("SJAR", 2457620.8,2457621.5,"R"),c("DUBF", 2457737.5, 2457738.8,"R"),c("LDJ",2457392.44315,2457394.44315,"V"))
editUser <- rbind(editUser,c("JM",2457466.9,2457468.6,"V"),c("PXR",2457573.0,2457574.5,"V"),c("MJB",2457634.0,2457635.3,"V"))
editUser <- rbind(editUser,c("HDHA",2457907,2457908.1,"V"),c("UJHA",2457925.2,2457926.5,"V"),c("SJAR",2457620.9,2457621.91,"I"),c("PALE",2457667.1,2457668.2,"V"))
editUser <- rbind(editUser,c("PXR",2457531.0,2457532.2,"B"),c("MJB",2457570.2, 2457571.2,"V"),c("DUBF",2457628.8,2457629.9,"B"))
editUser <- rbind(editUser,c("HBB", 2457489.3, 2457490.4,"V"),c("KTHC",2457956.0,2457957.2,"V"),c("MJB",2457624.2,2457625.3,"V"),c("MJB",2457660.2,2457661.6,"V"))
editUser <- rbind(editUser,c("SJAR", 2457620.9, 2457622.0,"V"),c("BJFB",2457925,2457928,"V"),c("LPB",2457348,2457349.5,"B"),c("LPB",2457385,2457386,"R"),
							c("LPB",2457348,2457349.5,"V"),c("SJAR",2457363,2457365,"R"),c("SJAR",2457363,2457365,"V"),c("SJAR",2458013.7,2458014.9,"R"),
							c("JM",2457982,2457983.4,"R"),c("JM",2457717,2457718.2,"R"),c("JM",2457468.4,2457469.5,"R"),c("SGEA",2457929.3, 2457930.4,"V"),
							c("JM",2457886.2,2457887.4,"B"),c("JM",2457995.2,2457996.3,"R"),c("OJJ",2458037.2, 2458038.2,"I"),
                            c("JSJA",2457899,2457900,"B"),c("JM",2458070.1,2458071.1,"B"),c("JM", 2458054.1, 2458055.1,"B"),c("OJJ",2458042.1,2458043.2,"I"),
                            c("DUBF",2457680	,2457790,"R"),c("JM",2458000,2458003,"B"),c("JM",2458048.1,2458049.2,"V"),c("JM",2457731.5,2457732.6,"B"),
                            c("AMID",2458073.8,2458074.8,"V"),c("AMID",2457937.0,2457938.0,"V"),c("ASASSN",2457689.2,2457690.3,"V"),c("JM",2458317.69,2458317.77,"R"))

############################################################
use.static.biases <- TRUE
biasObserver <- data.frame(obsCode=character(),band=character(),bias=numeric(),stringsAsFactors=FALSE)
# list each observer for each band only once
biasObserver <- rbind(biasObserver, list(obsCode="OAR",band="B",bias=-0.0125),
									list("JM","B",-0.008),
                                    list("JM","V",-0.018),
									list("JM","R",0.0),
									list("JM","I",0.017),
									list("LDJ","R",0.0),
									list("LDJ","I",0.003),
									list("LDJ","B",0.007),
                                    list("LDJ","V",-0.004),
                                    list("LPB","R",0.07),
                                    list("GKA","R",0.015),
                                    list("GKA","I",0.014),
                                    list("GKA","V",0.005),
                                    list("HDHA","V",0.01),
                                    list("HDHA","I",-0.015),
                                    list("SDB","V",-0.015),
                                    list("PXR","V",-0.010),
                                    list("OAR","V",0.011),
                                    list("OAR","B",-0.007),
                                    list("OAR","I",0.018),
                                    list("HJW","V",0.016),
                                    list("HJW","B",0.036),
                                    list("DUBF","B",-0.003),
                                    list("DUBF","V",+0.0034),
                                    list("DUBF","R",-0.029),
                                    list("DKS","B",-0.0015),
                                    list("MJB","I",0.01),
                                    list("MJB","V",0.007),
                                    list("MJB","B",0.002),
                                    list("OJJ","V",-0.006),
                                    list("OJJ","I",-0.039),
                                    list("LBG","I",0.023),
                                    list("LBG","R",0.009),
                                    list("LBG","B",0.03),
                                    list("JSJA","V",-0.004),
                                    list("SDB","V",0.005),
                                    list("SGEA","V",0.027),
                                    list("LPB","R",-0.013),
                                    list("CPP","B",0.027),
                                    list("LPAC","B",0.03),
                                    list("BMAK","B",-0.01),
                                    list("LWHA","B",0.0155),
                                    list("HBB","V",0.019),
                                    list("HBB","B",0.048),
                                    list("ASASSN","V",0.04))
