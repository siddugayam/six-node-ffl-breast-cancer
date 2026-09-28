import asyncio, re
from playwright.async_api import async_playwright
SP="/path/to/scratch"
URL="https://kmplot.com/analysis/index.php?p=service&cancer=breast"
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(); pg=await b.new_page()
        bodies=[]
        pg.on("request", lambda r: bodies.append((r.url,r.post_data)) if (r.method=="POST" and "kmplot" in r.url) else None)
        await pg.goto(URL, wait_until="networkidle", timeout=120000)
        # survival options
        surv=await pg.evaluate("Array.from(document.querySelectorAll('#surv option')).map(o=>o.value+' | '+o.textContent.trim())")
        print("SURV OPTIONS:", surv)
        for sel in ['er_status_ihc','er_status_array','pr_status_ihc','her2_status_ihc_fish','her2_status_array','intrinsic_subtype','molecular_subtype_stgallen','grade','lymph_node_status','pietenpol_subtype','tp53_status','histology']:
            try:
                o=await pg.evaluate(f"Array.from(document.querySelectorAll('#{sel} option')).map(o=>o.value+' | '+o.textContent.trim())")
                if o: print("SELECT",sel,":",o)
            except Exception: pass
        await pg.fill("#affyid","PTEN")
        await pg.keyboard.type(" "); await pg.keyboard.press("Backspace")
        await pg.wait_for_timeout(3000)
        try:
            await pg.click("#affyid_autocomplete li:first-child", timeout=8000); print("clicked ac")
        except Exception as e: print("ac fail",e)
        await pg.wait_for_timeout(1500)
        await pg.check("#termsOfUseAccepted")
        await pg.evaluate("var e=document.getElementById('display_results_in_new_window'); if(e){e.checked=false;} document.getElementById('myForm').target='_self';")
        await pg.evaluate("document.querySelectorAll('#export_data_as_txt').forEach(e=>e.checked=true)")
        await pg.click("#start_button")
        await pg.wait_for_load_state("networkidle", timeout=240000)
        await pg.wait_for_timeout(3000)
        out=await pg.content()
        open(SP+"/km_mrna_result.html","w").write(out)
        print("len",len(out))
        for u,d in bodies:
            if d: open(SP+"/km_mrna_postbody.txt","w").write(d); print("captured body from",u)
        await b.close()
asyncio.run(main())
