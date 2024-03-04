from randomizerMain import BigRandomizer, RandLauncher
## Domo how to launch the randomization
templ = {
    "useLinkList": True,
    "inputFn": "./sampleOriginalLinks.csv",
    "outputFn": "./myRandomizedLinks.csv",
    "q": 0.7,
    "runNaiveRandomization": False,

    "coolRelated":{
    
    }
}
R = RandLauncher(templateName=templ)

#Or alternatively
#R = RandLauncher(templateName="templateLinks.json")
