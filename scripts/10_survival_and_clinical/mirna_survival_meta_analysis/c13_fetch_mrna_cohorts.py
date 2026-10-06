#!/usr/bin/env python3
# Downloads the GEO series matrices of the mRNA cohorts into cache/v6/cohorts/matrix/, skipping files already
# present.
import os,sys,time,urllib.request,urllib.error
D="/path/to/revision/cache/v6/cohorts/matrix"
os.makedirs(D,exist_ok=True)
def dl(url,path):
    if os.path.exists(path) and os.path.getsize(path)>1000:
        print("cached",os.path.basename(path)); return True
    for i in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"mirna-ffl/1.0"})
            with urllib.request.urlopen(req,timeout=900) as r,open(path,"wb") as f:
                while True:
                    b=r.read(1<<20)
                    if not b: break
                    f.write(b)
            print("OK",url,os.path.getsize(path),flush=True); return True
        except urllib.error.HTTPError as e:
            print("HTTP",e.code,url,flush=True); return False
        except Exception as e:
            print("retry",e,flush=True); time.sleep(5+5*i)
    return False
for gse,gpl in [("GSE1456","GPL96"),("GSE1456","GPL97"),("GSE7390",None),
                ("GSE11121",None),("GSE2034",None),("GSE110651",None),
                ("GSE118782",None),("GSE44281",None)]:
    stub=gse[:-3]+"nnn"
    fn="%s_series_matrix.txt.gz"%gse if gpl is None else "%s-%s_series_matrix.txt.gz"%(gse,gpl)
    dl("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s"%(stub,gse,fn), os.path.join(D,fn))
