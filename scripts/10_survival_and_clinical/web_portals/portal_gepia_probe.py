import asyncio
from playwright.async_api import async_playwright
SP="/path/to/scratch"
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await (await b.new_context(viewport={"width":1500,"height":1100})).new_page()
        reqs=[]
        pg.on("request", lambda r: reqs.append((r.method,r.url,r.post_data)))
        await pg.goto("http://gepia2.cancer-pku.cn/#survival", wait_until="networkidle", timeout=180000)
        await pg.wait_for_timeout(4000)
        html=await pg.content(); open(SP+"/gepia_surv_dom.html","w").write(html)
        ctrl=await pg.evaluate("""Array.from(document.querySelectorAll('input,select,button,textarea')).map(e=>{
            const r=e.getBoundingClientRect(); return {t:e.tagName,ty:e.type,id:e.id,name:e.name,val:e.value,vis:(r.width>0&&r.height>0)};})""")
        for c in ctrl:
            if c['vis']: print(c)
        print("=== requests ===")
        for m,u,d in reqs:
            if 'gepia' in u and ('php' in u or 'api' in u or m=='POST'): print(m,u,'|',(d or '')[:200])
        await b.close()
asyncio.run(main())
