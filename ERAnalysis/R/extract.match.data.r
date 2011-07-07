extract.match.data <- function ()
{
# Extracts the relevant data from recent surveys to do matching of double observer data.
# There are no arguments but it uses the dataframe ERSruveyData and Height which
# contains the platform heights at the observation locations each year.
# Returns a list of dataframes that are divided based on survey date.
#
# added ret, angle, dist to sighting, dist to Azmuth
#
    DistAzmuth <- function (sangle, sdist){return((sin(pi * (sangle - 241)/180) * sdist))}
    data(ERSurveyData,package="ERAnalysis",envir=environment())
    ERSurveyData$LOCATION[ERSurveyData$LOCATION=="NP"]="N"
    Height=NULL
    data(Height,package="ERAnalysis",envir=environment())
# Extract primary observer sightings in which travel direction is not North
# Primary observers are at South location except in years 2000,2001
# This could be used to create Primary dataframe that is in package.  Currently
# Primary and SecondarySightings are constructed from ERAbund files.
    xx = ERSurveyData[ERSurveyData$EFLAG == 3 & ERSurveyData$TRAVELDIR !=
        "N", ]
    xx = rbind(xx[xx$Start.year %in% c(2000, 2001) & (xx$EXPERIMENT ==
        1 | (xx$EXPERIMENT == 2 & xx$LOCATION == "N")), ], xx[!(xx$Start.year %in%
        c(2000, 2001)) & (xx$EXPERIMENT == 1 | (xx$EXPERIMENT ==
        2 & xx$LOCATION == "S")), ])
    Primaryx = xx
    Primaryx$station = "P"
# Extract secondary observer sightings in which travel direction is not North
# Secondary observers are at North location except in years 2000,2001
    xx = ERSurveyData[ERSurveyData$EFLAG == 3 & ERSurveyData$TRAVELDIR !=
        "N", ]
    xx = rbind(xx[xx$Start.year %in% c(2000, 2001) & (xx$EXPERIMENT ==
        2 & xx$LOCATION == "S"), ], xx[!(xx$Start.year %in% c(2000,
        2001)) & (xx$EXPERIMENT == 2 & xx$LOCATION != "S"), ])
    Secondaryx = xx
    Secondaryx$station = "S"
# Create Match.data dataframe with the experiment 2 records
    Match.data = rbind(Primaryx[Primaryx$EXPERIMENT == 2, ],
        Secondaryx[Secondaryx$EXPERIMENT == 2, ])
# get reticles,angles and times using South sightings if available otherwise use North
    reticles = Match.data$SRET
    reticles[floor(reticles * 1000) == 0] = Match.data$NRET[floor(reticles *
        1000) == 0]
    reticles[floor(reticles * 1000) == 0] = NA
    angles = Match.data$SANGLE
    angles[floor(angles) == 0] = Match.data$NANGLE[floor(angles) == 0]
    times = Match.data$STIME
    times[floor(times) == 0] = Match.data$NTIME[floor(times) ==  0]
# constuct key field to match into Height dataframe to get heights for
# distance calculation from reticles
    Match.data$key = paste(Match.data$Start.year,Match.data$LOCATION,sep = "")
    heights = merge(Match.data, Height, all.x = TRUE, by = "key")$Height
    distances = RetDistK(Height = heights, Reticles = reticles)
# convert radial distance to perpendicular distance
    Match.data$sdup = RetDistK(Height = heights, Reticles = (reticles-.05))
    Match.data$sddn = RetDistK(Height = heights, Reticles = (reticles+.05))
    Match.data$sdist = distances
    Match.data$sang2A = angles-241
    Match.data$stime = times
    Match.data$sret = reticles
    Match.data$distance = D241(angles, distances)
    Match.data$sd2A = DistAzmuth(angles, distances)
# Compute crossing time at the line perpendicular to the coast at the shore station (abeam)
    Match.data$t241 = T241H(times, angles, distances, Speed = 6)
# Compute day since 1 Dec
    Match.data$day = as.vector(Match.data$DATE - as.POSIXct(paste(Match.data$Start.year,
        "-12-01", sep = ""))+1)
# Compute watch factor variable
    Match.data$watch = 1
    Match.data$watch[Match.data$t241 >= 10.5] = 2
    Match.data$watch[Match.data$t241 >= 13.5] = 3
# Create dataframe excluding those with missing T241 or distance and modify field names
    t241=NULL
    distance=NULL
    Match.data = subset(Match.data, subset=(t241 & distance),
              select = c("station", "day",  "watch", "t241", "distance"
        , "PODSIZE","stime", "sret", "sdist", "sang2A", "sd2A","sdup","sddn"
        , "OBSERVER", "VISCODE", "WINDFORCE", "WINDDIR"
        , "Start.year"))
    names(Match.data)[names(Match.data) == "PODSIZE"] = "podsize"
   names(Match.data)[names(Match.data) == "OBSERVER"] = "observer"
   names(Match.data)[names(Match.data) == "VISCODE"] = "vis"
    names(Match.data)[names(Match.data) == "WINDFORCE"] = "beaufort"
    names(Match.data)[names(Match.data) == "WINDDIR"] = "wind.direction"
 # The following will spit the Match by DATE to do the matching.
    Match.list = split(Match.data, paste(Match.data$Start.year,
        Match.data$day, sep = "_"))
    rm(ERSurveyData)
    return(Match.list)
}





 
