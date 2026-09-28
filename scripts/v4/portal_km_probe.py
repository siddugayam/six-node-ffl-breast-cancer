import sys, json, re, asyncio
from playwright.async_api import async_playwright

URL="https://kmplot.com/analysis/index.php?p=service&cancer=breast_mirna"

async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch()
        pg=await b.new_page()
        posts=[]
        pg.on("request", lambda r: posts.append((r.method,r.url,r.post_data)) if r.method=="POST" else None)
        await pg.goto(URL, wait_until="networkidle", timeout=120000)
        await pg.fill("#affyid","hsa-miR-130a")
        await pg.keyboard.type(" ")
        await pg.keyboard.press("Backspace")
        await pg.wait_for_timeout(3000)
        html=await pg.content()
        open("/path/to/scratch/ac_dom.html","w").write(html)
        # try clicking autocomplete item
        try:
            await pg.click("#affyid_autocomplete li:first-child", timeout=8000)
            print("clicked autocomplete")
        except Exception as e:
            print("no autocomplete click:", e)
        await pg.wait_for_timeout(2000)
        print("affyid_hidden=", await pg.input_value("#affyid_hidden"))
        gl = await pg.evaluate("document.getElementById('genes_list') ? document.getElementById('genes_list').outerHTML : 'NO genes_list'")
        print("genes_list:", gl[:2000])
        # check terms
        await pg.check("#termsOfUseAccepted")
        await pg.evaluate("var e=document.getElementById('display_results_in_new_window'); if(e){e.checked=false;} var f=document.getElementById('myForm'); if(f){f.target='_self';}")
        await pg.evaluate("document.querySelectorAll('#export_data_as_txt').forEach(e=>e.checked=true)")
        # dataset METABRIC
        await pg.check("#METABRIC")
        await pg.click("#start_button")
        await pg.wait_for_load_state("networkidle", timeout=180000)
        await pg.wait_for_timeout(3000)
        out=await pg.content()
        open("/path/to/scratch/km_browser_result.html","w").write(out)
        print("RESULT URL:", pg.url, "len", len(out))
        for m,u,d in posts:
            print("POST",u)
            if d: open("/path/to/scratch/km_postbody.txt","w").write(d)
        await b.close()
asyncio.run(main())
