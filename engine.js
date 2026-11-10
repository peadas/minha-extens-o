/* ===== motor de dados compartilhado (ficha + sobreposição) =====
   Física de corpo rígido determinística: mesma "semente" => mesma jogada em qualquer computador.
   Não usa Math.hypot/sin/cos/acos dentro da simulação (podem variar de navegador para navegador). */
(function(root){
"use strict";

/* ---------- semente ---------- */
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296}}

/* ---------- física ---------- */
const PH=(()=>{
const P={G:2800,REST_F:.52,REST_W:.5,REST_D:.45,MU:.55,THR:120,IT:6,PC:.4,LIN_D:.06,ANG_D:.15,FLOOR_LD:2.4,FLOOR_AD:3,SV:30,SW:1.6,TILT_COS:0.9855847669095608 /* cos(0,17 rad) */};
const add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],mul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len=a=>Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
const qmul=(a,b)=>[a[3]*b[0]+a[0]*b[3]+a[1]*b[2]-a[2]*b[1],a[3]*b[1]-a[0]*b[2]+a[1]*b[3]+a[2]*b[0],a[3]*b[2]+a[0]*b[1]-a[1]*b[0]+a[2]*b[3],a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2]];
const qrot=(q,v)=>{const u=[q[0],q[1],q[2]],t=mul(cross(u,v),2);return add(add(v,mul(t,q[3])),cross(u,t))};
const qinv=q=>[-q[0],-q[1],-q[2],q[3]];
const qnorm=q=>{const l=Math.sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);return[q[0]/l,q[1]/l,q[2]/l,q[3]/l]};
function randQuat(rnd){let a,b,c,d,l;do{a=rnd()*2-1;b=rnd()*2-1;c=rnd()*2-1;d=rnd()*2-1;l=a*a+b*b+c*c+d*d}while(l>1||l<.01);l=Math.sqrt(l);return[a/l,b/l,c/l,d/l]}

/* hull = {verts:[[x,y,z]], planes:[{n:[x,y,z], d, vi:[índices]}]} (em unidades locais) */
class Body{
  constructor(hull,s){
    this.verts=hull.verts.map(v=>mul(v,s));
    this.planes=hull.planes.map(p=>({n:p.n,d:p.d*s,vi:p.vi}));
    this.rad=Math.max(...this.verts.map(len));
    this.invI=1/(.25*this.rad*this.rad);
    this.p=[0,0,0];this.q=[0,0,0,1];this.v=[0,0,0];this.w=[0,0,0];this.touch=false;
  }
}
class World{
  constructor(b){this.setBounds(b);this.bodies=[];this.slowT=0;this.t=0}
  setBounds(b){
    this.pl=[[[0,1,0],0,true],[[1,0,0],b.minX,false],[[-1,0,0],-b.maxX,false],[[0,0,1],b.minZ,false],[[0,0,-1],-b.maxZ,false]];
  }
  step(dt){
    const B=this.bodies;
    for(const b of B){
      b.v[1]-=P.G*dt;
      b.v=mul(b.v,1-P.LIN_D*dt);b.w=mul(b.w,1-P.ANG_D*dt);
      b.p=add(b.p,mul(b.v,dt));
      const dq=qmul([b.w[0],b.w[1],b.w[2],0],b.q);
      b.q=qnorm([b.q[0]+.5*dt*dq[0],b.q[1]+.5*dt*dq[1],b.q[2]+.5*dt*dq[2],b.q[3]+.5*dt*dq[3]]);
      b.touch=false;
    }
    for(let it=0;it<P.IT;it++){
      for(const b of B)this.planes(b);
      for(let i=0;i<B.length;i++)for(let j=i+1;j<B.length;j++)this.pair(B[i],B[j]);
    }
    let slow=true;
    for(const b of B){
      if(b.touch){b.v=mul(b.v,1-P.FLOOR_LD*dt);b.w=mul(b.w,1-P.FLOOR_AD*dt)}
      if(len(b.v)>P.SV||len(b.w)>P.SW)slow=false;
    }
    this.slowT=slow?this.slowT+dt:0;this.t+=dt;
  }
  done(){return this.slowT>=.5}
  planes(b){
    for(const [n,off,fl] of this.pl){
      let maxPen=0;
      for(let i=0;i<b.verts.length;i++){
        const r=qrot(b.q,b.verts[i]),pen=off-dot(n,add(b.p,r));
        if(pen<=0)continue;
        if(pen>maxPen)maxPen=pen;
        if(fl)b.touch=true;
        let vel=add(b.v,cross(b.w,r));const vn=dot(vel,n);
        if(vn<0){
          const e=-vn>P.THR?(fl?P.REST_F:P.REST_W):0,rn=cross(r,n),j=-(1+e)*vn/(1+b.invI*dot(rn,rn));
          b.v=add(b.v,mul(n,j));b.w=add(b.w,mul(rn,j*b.invI));
          vel=add(b.v,cross(b.w,r));
          const vt=sub(vel,mul(n,dot(vel,n))),vl=len(vt);
          if(vl>1e-4){
            const t=mul(vt,1/vl),rt=cross(r,t),jt=Math.min(vl/(1+b.invI*dot(rt,rt)),P.MU*j);
            b.v=sub(b.v,mul(t,jt));b.w=sub(b.w,mul(rt,jt*b.invI));
          }
        }
      }
      if(maxPen>0)b.p=add(b.p,mul(n,maxPen*P.PC));
    }
  }
  pair(a,b){
    if(len(sub(b.p,a.p))>a.rad+b.rad)return;
    this.vin(a,b);this.vin(b,a);
  }
  /* vértices de a dentro do casco de b */
  vin(a,b){
    const qi=qinv(b.q);
    for(const vl of a.verts){
      const r=qrot(a.q,vl),x=add(a.p,r),xl=qrot(qi,sub(x,b.p));
      let pen=1e9,fi=-1;
      for(let k=0;k<b.planes.length;k++){const dd=b.planes[k].d-dot(b.planes[k].n,xl);if(dd<pen){pen=dd;fi=k}}
      if(pen<=0)continue;
      const n=qrot(b.q,b.planes[fi].n),rb=sub(x,b.p);
      const rel=sub(add(a.v,cross(a.w,r)),add(b.v,cross(b.w,rb))),vn=dot(rel,n);
      if(vn<0){
        const e=-vn>P.THR?P.REST_D:0,ra=cross(r,n),rbn=cross(rb,n);
        const j=-(1+e)*vn/(2+a.invI*dot(ra,ra)+b.invI*dot(rbn,rbn));
        a.v=add(a.v,mul(n,j));a.w=add(a.w,mul(ra,j*a.invI));
        b.v=sub(b.v,mul(n,j));b.w=sub(b.w,mul(rbn,j*b.invI));
        const r2=sub(add(a.v,cross(a.w,r)),add(b.v,cross(b.w,rb))),vt=sub(r2,mul(n,dot(r2,n))),vl2=len(vt);
        if(vl2>1e-4){
          const t=mul(vt,1/vl2),rat=cross(r,t),rbt=cross(rb,t),jt=Math.min(vl2/(2+a.invI*dot(rat,rat)+b.invI*dot(rbt,rbt)),P.MU*j);
          a.v=sub(a.v,mul(t,jt));a.w=sub(a.w,mul(rat,jt*a.invI));
          b.v=add(b.v,mul(t,jt));b.w=add(b.w,mul(rbt,jt*b.invI));
        }
      }
      const c=pen*.25;a.p=add(a.p,mul(n,c));b.p=sub(b.p,mul(n,c));
    }
  }
}
/* lança o dado a partir do meio da arena; i/n = posição na fila para não nascerem sobrepostos */
function toss(b,ar,rnd,i,n){
  const cx=(ar.minX+ar.maxX)/2,cz=(ar.minZ+ar.maxZ)/2,sp=Math.min(b.rad*2.2,(ar.maxX-ar.minX-4*b.rad)/Math.max(1,n-1));
  /* caem de cima, juntos, quase um em cima do outro: batem entre si e se espalham */
  b.p=[cx+(i-(n-1)/2)*sp+(rnd()-.5)*b.rad*1.2,300+rnd()*160+i*(b.rad*2.6),cz+(rnd()-.5)*b.rad*2.2];
  /* direção horizontal sorteada sem trigonometria (resultado idêntico em qualquer navegador) */
  let x,z,l;do{x=rnd()*2-1;z=rnd()*2-1;l=x*x+z*z}while(l>1||l<.01);l=Math.sqrt(l);
  const s=800+rnd()*900;
  b.v=[x/l*s,-100-rnd()*200,z/l*s];
  const ax=[rnd()*2-1,rnd()*2-1,rnd()*2-1],ln=len(ax)||1,m=22+rnd()*26;
  b.w=mul(ax,m/ln);b.q=randQuat(rnd);
}
/* dado "torto" (apoiado na parede, em outro dado…): pega e joga de novo */
function kick(b,rnd){
  b.v=[(rnd()-.5)*300,320+rnd()*220,(rnd()-.5)*300];
  const ax=[rnd()*2-1,rnd()*2-1,rnd()*2-1],l=len(ax)||1;b.w=mul(ax,(14+rnd()*10)/l);
}
/* lê o resultado: face voltada para cima (d4: vértice voltado para cima); bad = está torto */
function readout(b,d4){
  let bi=0,best=d4?1e9:-1e9;
  b.planes.forEach((p,i)=>{const y=qrot(b.q,p.n)[1];if(d4?y<best:y>best){best=y;bi=i}});
  const bad=Math.abs(best)<P.TILT_COS;
  if(!d4)return{idx:bi,bad};
  const used=b.planes[bi].vi;let o=0;
  for(let i=0;i<b.verts.length;i++)if(!used.includes(i)){o=i;break}
  return{idx:o,bad};
}
/* orientação que faz o resultado "want" (índice 0..f-1) ficar virado para cima */
function targetQuat(b,want,d4){
  let k=want,dir=[0,1,0];
  if(d4){k=b.planes.findIndex(p=>!p.vi.includes(want));if(k<0)k=0;dir=[0,-1,0]}
  const cur=qrot(b.q,b.planes[k].n),d=dot(cur,dir);let r;
  if(d<-.999999){let ax=cross(cur,[1,0,0]);if(len(ax)<1e-6)ax=cross(cur,[0,0,1]);ax=mul(ax,1/len(ax));r=[ax[0],ax[1],ax[2],0]}
  else{const c=cross(cur,dir);r=qnorm([c[0],c[1],c[2],1+d])}
  return qnorm(qmul(r,b.q));
}
return{P,Body,World,toss,kick,readout,targetQuat,qrot,randQuat};
})();

/* ---------- geometria (precisa do three.js, só na hora de usar) ---------- */
const GEO={
 4:()=>new THREE.TetrahedronGeometry(1.15),
 6:()=>new THREE.BoxGeometry(1.35,1.35,1.35),
 8:()=>new THREE.OctahedronGeometry(1.1),
 10:()=>{const p=[],T=[0,1.1,0],B=[0,-1.1,0],U=[],L=[],r=.95,h=.116;
  for(let i=0;i<5;i++){const a=i*Math.PI*2/5,b=a+Math.PI/5;U.push([Math.cos(a)*r,h,Math.sin(a)*r]);L.push([Math.cos(b)*r,-h,Math.sin(b)*r])}
  for(let i=0;i<5;i++){const j=(i+1)%5;p.push(...T,...U[i],...L[i],...T,...L[i],...U[j],...B,...L[i],...U[j],...B,...U[j],...L[j])}
  const g=new THREE.BufferGeometry();g.setAttribute("position",new THREE.Float32BufferAttribute(p,3));g.computeVertexNormals();return g},
 12:()=>new THREE.DodecahedronGeometry(1.05),
 20:()=>new THREE.IcosahedronGeometry(1.1)
};
function faces(g){
  const p=g.attributes.position,F=[],key=v=>v.toArray().map(x=>x.toFixed(3)).join();
  for(let i=0;i<p.count;i+=3){
    let v=[0,1,2].map(k=>new THREE.Vector3().fromBufferAttribute(p,i+k));
    const n=new THREE.Vector3().subVectors(v[1],v[0]).cross(new THREE.Vector3().subVectors(v[2],v[0])).normalize();
    if(n.dot(v[0].clone().add(v[1]).add(v[2]))<0){v=[v[0],v[2],v[1]];n.negate()}
    let f=F.find(f=>f.n.dot(n)>.995);
    if(!f){f={n,t:[],u:{}};F.push(f)}
    f.t.push(v);v.forEach(w=>f.u[key(w)]=w);
  }
  F.forEach(f=>{const u=Object.values(f.u);f.c=u.reduce((a,w)=>a.add(w),new THREE.Vector3()).divideScalar(u.length);f.R=u.reduce((a,w)=>a+w.distanceTo(f.c),0)/u.length});
  return F;
}
const cache={};
/* analisa o poliedro de um dado: faces, vértices e o "casco" usado pela física */
const ENG={PH,GEO,faces,mulberry32};
ENG.analyze=function(f){
  if(cache[f])return cache[f];
  const g=GEO[f](),Fs=faces(g.index?g.toNonIndexed():g),vm={},V=[];
  Fs.forEach(F=>{F.vi=Object.keys(F.u).map(k=>{if(!(k in vm)){vm[k]=V.length;V.push(F.u[k])}return vm[k]})});
  const mx=Math.max(...V.map(v=>v.length()));
  const hull={verts:V.map(v=>v.toArray()),planes:Fs.map(F=>({n:F.n.toArray(),d:F.n.dot(F.c),vi:F.vi}))};
  return cache[f]={Fs,V,hull,mx};
};

/* ---------- arena fixa (igual para todo mundo; a tela só muda o zoom) ---------- */
ENG.ARENA={minX:-540,maxX:540,minZ:-300,maxZ:300};
ENG.RAD=38;
ENG.VIEW={w:1140,h:660};
ENG.DT=1/120;
ENG.FACES=[4,6,8,10,12,20];

/* ---------- simulação passo a passo ---------- */
class Sim{
  constructor(list,seed){
    this.list=list;this.rnd=mulberry32(seed);this.ar=ENG.ARENA;
    this.W=new PH.World(this.ar);this.bodies=[];this.kicks=0;this.finished=false;this.results=null;this.flat=true;
    list.forEach(([n,f],i)=>{
      const a=ENG.analyze(f),b=new PH.Body(a.hull,ENG.RAD/a.mx);
      b.f=f;PH.toss(b,this.ar,this.rnd,i,list.length);this.W.bodies.push(b);this.bodies.push(b);
    });
  }
  stepOnce(){
    if(this.finished)return;
    const W=this.W,ar=this.ar;
    W.step(ENG.DT);
    for(const b of this.bodies){
      if(!isFinite(b.p[0]+b.p[1]+b.p[2]))throw new Error("simulação instável");
      b.p[0]=Math.min(ar.maxX,Math.max(ar.minX,b.p[0]));b.p[2]=Math.min(ar.maxZ,Math.max(ar.minZ,b.p[2]));b.p[1]=Math.min(1500,b.p[1]);
    }
    if(W.done()||W.t>8){
      const rd=this.bodies.map(b=>PH.readout(b,b.f===4));
      if(rd.some(r=>r.bad)&&this.kicks<3&&W.t<8){
        this.kicks++;rd.forEach((r,i)=>{if(r.bad)PH.kick(this.bodies[i],this.rnd)});W.slowT=0;
      }else{
        this.finished=true;this.flat=!rd.some(r=>r.bad);
        this.results=this.list.map(([n,f],i)=>[n,f,rd[i].idx+1]);
      }
    }
  }
  /* roda tudo de uma vez (sem desenhar): é assim que a ficha descobre o resultado */
  runAll(){let n=0;while(!this.finished&&n++<4000)this.stepOnce();if(!this.finished)throw new Error("simulação não terminou");return this.results}
}
ENG.Sim=Sim;

/* ---------- validação de dados vindos de outros jogadores ---------- */
const str=(s,n)=>String(s==null?"":s).replace(/[\u0000-\u001f\u007f]/g," ").slice(0,n);
const KIND=["","pv","rel"],CLS=["x","pv","rl","ok","bad"];
ENG.clean=function(p){
  if(!p||typeof p!=="object"||!Array.isArray(p.list))return null;
  const list=[],kinds=[],res=[];
  for(const x of p.list.slice(0,8)){
    if(!Array.isArray(x))continue;
    const f=+x[1];if(!ENG.FACES.includes(f))continue;
    const i=list.length;
    list.push([str(x[0],24),f]);
    kinds.push(KIND.includes(p.kinds&&p.kinds[i])?p.kinds[i]:"");
    const r=Math.floor(+(p.res&&p.res[i]));res.push(r>=1&&r<=f?r:1);
  }
  if(!list.length)return null;
  const cls=s=>str(s,20).split(/\s+/).filter(t=>CLS.includes(t)).join(" ");
  const rows=(Array.isArray(p.rows)?p.rows:[]).slice(0,14).map(r=>({l:str(r&&r.l,60),v:str(r&&r.v,12),c:cls(r&&r.c)}));
  const crit=p.crit&&typeof p.crit==="object"?{t:str(p.crit.t,50),c:cls(p.crit.c)}:null;
  return{id:str(p.id,32),seed:(+p.seed)>>>0,who:str(p.who,40)||"Jogador",color:/^#[0-9a-f]{6}$/i.test(p.color)?p.color:"#ecc962",
    title:str(p.title,40),list,kinds,res,anim:p.anim!==false,rows,crit,note:p.note?str(p.note,120):"",text:str(p.text,160)};
};
ENG.b64=o=>btoa(unescape(encodeURIComponent(JSON.stringify(o)))).replace(/\+/g,"-").replace(/\//g,"_").replace(/=+$/,"");
ENG.unb64=s=>JSON.parse(decodeURIComponent(escape(atob(String(s).replace(/-/g,"+").replace(/_/g,"/")))));
ENG.CH={roll:"ordem.ficha/roll",feed:"ordem.ficha/feed",idle:"ordem.ficha/idle"};

root.ENG=ENG;
if(typeof module!=="undefined"&&module.exports)module.exports=ENG;
})(typeof window!=="undefined"?window:globalThis);
