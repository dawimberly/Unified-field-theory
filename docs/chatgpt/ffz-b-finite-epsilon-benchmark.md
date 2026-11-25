# FFff(B) unsubtracted finite-epsilon benchmark

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records FFff(B) unsubtracted finite-epsilon benchmark. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues the [homogeneous overlap audit](ffz-b-homogeneous-tail-audit.md). This evaluates the full projective kernel in its convergent domain without any Laurent coefficients as inputs.

**ASSUMED** Set the effective perimeter to one and use $Q=(1-P)q$. The scalar and fermion physical lengths of the first mixed pair are $P/a$ and $P/(1-a)$. Include $P^{\epsilon -1} Q^{\epsilon -1}$, the two fraction weights $[a(1-a)]^{-\epsilon }(1-a)$ and $[b(1-b)]^{-\epsilon }(1-b)$, the bridge Jacobian 1−P−Q, and $d(d+2)/4$. The curvature term is $W^{-3\epsilon }(WA+kB)$, where W is the physical perimeter. The table excludes the outer radial gamma factor and curvature normalization.

**DERIVED** For the bridge division the moments have degree at most two and the weighted bridge sum degree at most three. Simpson integration is exact for these polynomials. The remaining quadrature removes the explicit P and q powers and clusters points near the possibly long constituent. At $\epsilon >1/2$ the isolated long scalar endpoint has integrable power $-2+2\epsilon$. Two anchor choices test the implementation in this domain.

| $\epsilon$ | Points per quadrature dimension | Anchor 0 | Anchor 2 |
| --- | ---: | ---: | ---: |
| **DERIVED** 0.7 | 40 | 1.0032901597 | 1.0032933363 |
| **DERIVED** 0.7 | 64 | 1.0033353814 | 1.0033362220 |
| **DERIVED** 0.8 | 40 | 0.5846625209 | 0.5846720891 |
| **DERIVED** 0.8 | 64 | 0.5849868907 | 0.5849900019 |

**DERIVED** These are numerical reference values. The order drift exceeds the anchor spread and limits the useful precision. The two anchor values converge toward one another. Continuation through the nonintegrable endpoints is required to test the double pole. This check does not resolve the mismatch.

## Reproduction

**DERIVED** Node.js, standard library.

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
function evaluate(n,eps,root=0){
 const d=4-2*eps,k=-1+3*eps, g=nodes(n);
 const pairs=g.flatMap(([x,w])=>{
  const q=.5*x**(1/(2*eps-1)), jac=q/((2*eps-1)*x);
  return [q,1-q].map(a=>[a,w*jac*(a*(1-a))**(-eps)*(1-a)]);
 });
 const lengths=g.map(([x,w])=>[x**(1/eps),w/eps]);
 let acc=0;
 for(const [P,wp] of lengths) for(const [q,wq] of lengths){
  const Q=(1-P)*q,r=(1-P)*(1-q);
  const m0=geometry(P,Q,0,d,root),m1=geometry(P,Q,.5,d,root),m2=geometry(P,Q,1,d,root);
  const [V,M,N,B]=m0.map((z,i)=>(z+4*m1[i]+m2[i])/6);
  let inner=0;
  for(const [a,wa] of pairs) for(const [b,wb] of pairs){
   const s=P/a,f=P/(1-a),t=Q/b,h=Q/(1-b),W=r+s+f+t+h;
   const A=-eps*V-d*(4-d)*(s+t)/3;
   const profile=2/3*(B+(s+f)*M+(t+h)*N+d*(s*s+t*t)+(d-1)*(f*f+h*h));
   inner+=wa*wb*W**(-3*eps)*(W*A+k*profile);
  }
  acc+=wp*wq*(1-P)**eps*r*inner;
 }
 return d*(d+2)/4*acc;
}
for(const e of [.7,.8]) for(const n of [12,24,40,64])
 console.log(e,n,evaluate(n,e),evaluate(n,e,2));
```

- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005). — eq./page UNVERIFIED.
- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
