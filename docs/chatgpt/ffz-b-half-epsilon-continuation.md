# FFff(B) continuation across epsilon = 1/2

**Conversation date range** 2026-11-22

**DERIVED** Subtracting the common long-endpoint model gives a finite $4\pi$ add-back at $\epsilon=1/2$. **DERIVED** Pointwise evaluation at that value misses the add-back. **UNVERIFIED** Quadrature drift and overlapping endpoints leave the double-pole discrepancy open.

## Research record

**ASSUMED** Continues the [unsubtracted benchmark](ffz-b-finite-epsilon-benchmark.md). This removes the first apparent singularity in continuing the projective integral toward $\epsilon =0$. It identifies a finite contribution that would be lost by setting $\epsilon =1/2$ directly.

**DERIVED** For the scalar-long endpoint $a\to 0$ at fixed P, the largest physical length is $s=P/a$. The scalar curvature-mass term in WA and the scalar variance term in kB contribute jointly a coefficient $2d(2\epsilon -1)/3$ times $s^{2}$, since $k=-1+3\epsilon$ and $d=4-2\epsilon$. The leading integrand in a is therefore $2d(2\epsilon -1)P^{2-3\epsilon }a^{2\epsilon -2}/3$, apart from the partner measure. Subtract this leading term for $a<1/2$ and for the corresponding endpoint in b. Its individual integral contains $1/(2\epsilon -1)$, which cancels the explicit zero in the numerator. The two add-backs after spin normalization are exactly

$$
2\frac{2d}{3}\frac{d(d+2)}{4}2^{1-2\epsilon}
\frac{\Gamma(2-2\epsilon)\Gamma(\epsilon)}{\Gamma(4-\epsilon)}
B(1-\epsilon,2-\epsilon).
$$

**DERIVED** At $\epsilon =1/2$ this gives **$4\pi$** in the projective normalization. It is an add-back, not the complete value of the integral. In particular, a raw pointwise evaluation at $\epsilon =1/2$ drops this finite contribution. The ratio singularity is removable in the continued integral.

**UNVERIFIED** The following quadrature results are exploratory. The subtraction smooths passage across $\epsilon =1/2$. At $\epsilon =0.7$ the continued representation and the unsubtracted one are algebraically the same, but the present low-order quadrature of the subtracted representation still has substantial drift and does not give a precise numerical match. Near $\epsilon =0$ the short and overlapping long logarithms still need explicit treatment. No conclusion about the double-pole discrepancy follows from this check.

| $\epsilon$ | 12-point rule | 24-point rule |
| --- | ---: | ---: |
| **DERIVED** 0.49 | 2.56034 | 2.57408 |
| **DERIVED** 0.50 | 2.47064 | 2.48397 |
| **DERIVED** 0.51 | 2.38108 | 2.39379 |
| **DERIVED** 0.70 | 1.28358 | 1.12473 |

## Reproduction

```javascript
function nodes(n){
 const out=[];
 for(let j=1;j<=n;j++){
  let x=Math.cos(Math.PI*(j-.25)/(n+.5)),d;
  for(let i=0;i<8;i++){
   let p=1,q=x;
   for(let k=2;k<=n;k++) [p,q]=[q,((2*k-1)*x*q-(k-1)*p)/k];
   d=n*(p-x*q)/(1-x*x); x-=q/d;
  }
  out.push([(1+x)/2,1/((1-x*x)*d*d)]);
 }
 return out;
}
function geometry(P,Q,v,d,root){
 const r=1-P-Q, l=[P,r*v,Q,r*(1-v)];
 const t=[0,l[0],l[0]+l[1],l[0]+l[1]+l[2]];
 const b=t.map((_,i)=>l.map((_,j)=>2*((j<i)-(j<root)-t[i]+t[root])));
 const mom=(i,j)=>{
  const b1=b[i],b2=b[j];
  const f=(a,c)=>-(b1[a]*b2[c]+b1[c]*b2[a])/2;
  const ds=Math.abs(t[i]-t[j]), di=Math.abs(t[i]-t[root]),dj=Math.abs(t[j]-t[root]);
  const C=di*(1-di)+dj*(1-dj)-ds*(1-ds);
  return d*C+(d*(f(0,1)+f(2,3)+f(0,3)+f(1,2))+(2-d)*(f(0,2)+f(1,3)))/(d+2);
 };
 const m=t.map((_,i)=>mom(i,i));
 const e=m.map((_,i)=>m[i]+m[(i+1)%4]+mom(i,(i+1)%4));
 return [m.reduce((a,b)=>a+b),e[0],e[2],l[1]*e[1]+l[3]*e[3]+(d-1)*(l[1]**2+l[3]**2)];
}

function gamma(x){
 const c=[676.5203681218851,-1259.1392167224028,771.3234287776531,-176.6150291621406,12.507343278686905,-.13857109526572012,9.984369578019572e-6,1.5056327351493116e-7];
 x--; let y=.9999999999998099; c.forEach((v,i)=>y+=v/(x+i+1));
 const t=x+7.5;return Math.sqrt(2*Math.PI)*t**(x+.5)*Math.exp(-t)*y;
}
function restored(e){
 const d=4-2*e;
 return 2*(2*d/3)*d*(d+2)/4*.5**(2*e-1)
 *gamma(2-2*e)*gamma(e)/gamma(4-e)
 *gamma(1-e)*gamma(2-e)/gamma(3-2*e);
}
function run(n,eps){
 const d=4-2*eps,k=-1+3*eps, g=nodes(n);
 const pairs=g.flatMap(([x,w])=>[x*x/2,1-x*x/2].map((a,i)=>[a,w*x*(a*(1-a))**(-eps)*(1-a),i,w*x]));
 const lengths=g.map(([x,w])=>[x**(1/eps),w/eps]);
 let acc=0;
 for(const [P,wp] of lengths) for(const [q,wq] of lengths){
  const Q=(1-P)*q,r=(1-P)*(1-q);
  const m0=geometry(P,Q,0,d,0),m1=geometry(P,Q,.5,d,0),m2=geometry(P,Q,1,d,0);
  const [V,M,N,B]=m0.map((z,i)=>(z+4*m1[i]+m2[i])/6);
  let inner=0;
  for(const [a,wa,ia,ja] of pairs) for(const [b,wb,ib,jb] of pairs){
   const s=P/a,f=P/(1-a),t=Q/b,h=Q/(1-b),W=r+s+f+t+h;
   const c=2*d*(2*eps-1)/3;
   const smooth=-eps*W*V-d*(4-d)*((s+t)*(r+f+h)+2*s*t)/3
       +k*2/3*(B+(s+f)*M+(t+h)*N+(d-1)*(f*f+h*h));
   const w=W**(-3*eps);
   inner+=wa*wb*w*(smooth+c*(s*s+t*t))
     -c*((ia===0?ja*wb*P**(2-3*eps)*a**(2*eps-2):0)
        +(ib===0?jb*wa*Q**(2-3*eps)*b**(2*eps-2):0));
  }
  acc+=wp*wq*(1-P)**eps*r*inner;
 }
 return d*(d+2)/4*acc+restored(eps);
}
if(Math.abs(restored(.5)-4*Math.PI)>1e-12) throw Error("normalization");
for(const e of [.49,.5,.51,.7]) for(const n of [12,24])
 console.log(e,n,run(n,e));
```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
