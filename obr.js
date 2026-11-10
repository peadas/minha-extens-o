/* Carrega o SDK oficial do Owlbear Rodeo direto de uma CDN (sem precisar de npm).
   Tenta duas CDNs, porque se uma estiver fora do ar a outra segura. */
const URLS=[
  "https://cdn.jsdelivr.net/npm/@owlbear-rodeo/sdk@3.1.0/+esm",
  "https://esm.sh/@owlbear-rodeo/sdk@3.1.0"
];
const pick=m=>{for(const c of [m&&m.default,m&&m.default&&m.default.default,m&&m.OBR,m])if(c&&typeof c.onReady==="function")return c;return null};
export async function loadOBR(){
  for(const u of URLS){
    try{const O=pick(await import(u));if(O)return O}
    catch(e){console.warn("[ficha] SDK do Owlbear não carregou de",u,e)}
  }
  return null;
}
/* true quando a página está mesmo dentro do Owlbear e o SDK já conversa com ele */
export function ready(OBR){
  return new Promise(res=>{
    if(!OBR||!OBR.isAvailable){res(false);return}
    if(OBR.isReady){res(true);return}
    OBR.onReady(()=>res(true));
  });
}
