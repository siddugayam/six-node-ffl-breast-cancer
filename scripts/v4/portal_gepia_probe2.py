import asyncio, json
from playwright.async_api import async_playwright
SP="/path/to/scratch"
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); ctx=await b.new_context(viewport={"width":1500,"height":1100})
        pg=await ctx.new_page()
        reqs=[]
        pg.on("request", lambda r: reqs.append((r.method,r.url,r.post_data)))
        await pg.goto("http://gepia2.cancer-pku.cn/#survival", wait_until="networkidle", timeout=180000)
        await pg.wait_for_timeout(4000)
        await pg.fill("#survival_signature","PTEN")
        # choose BRCA in the multi-select
        opts=await pg.evaluate("Array.from(document.querySelectorAll('#survival_tcgat option')).map(o=>o.value).slice(0,60)")
        print("dataset options sample:",opts[:12])
        await pg.select_option("#survival_tcgat", ["BRCA"])
        await pg.click("#survival_plot")
        await pg.wait_for_timeout(20000)
        html=await pg.content(); open(SP+"/gepia_surv_result.html","w").write(html)
        print("=== gepia reqs after submit ===")
        for m,u,d in reqs:
            if 'gepia' in u and ('PHP' in u or m=='POST'): print(m,u,'|',(d or '')[:300])
        # find result frame content
        fr=await pg.evaluate("document.getElementById('survival_output')?document.getElementById('survival_output').innerHTML.slice(0,1500):'no survival_output'")
        print("OUT:",fr[:1200])
        await b.close()
asyncio.run(main())
