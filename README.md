# JART-TUD VCM memristor model

## Memristor Current Equation
### Simplified memristor current ($I_{\rm M}$) expression

```math
 \begin{aligned}[b]
  I_{\rm M}\left(N_{\rm d},V_{\rm M}, d_{\rm r}, d_{\rm l}\right) = &\textcolor{green}{{p}_{\rm1}(V_{\rm M}, d_{\rm r}, d_{\rm l})}\left(\textcolor{green}{{p}_{\rm2}(V_{\rm M}, d_{\rm r}, d_{\rm l})}\left( e^{\frac{\left(\ln\frac{N_{\rm d}}{N_{\rm d,L}}-\textcolor{green}{{p}_{\rm3}(V_{\rm M}, d_{\rm r}, d_{\rm l})}\right)}{\textcolor{green}{{p}_{\rm4}(V_{\rm M}, d_{\rm r}, d_{\rm l})}}}-1\right)+\left(\ln\frac{N_{\rm d}}{N_{\rm d,L}}-\textcolor{green}{{p}_{\rm3}(V_{\rm M}, d_{\rm r}, d_{\rm l})}\right)\right) \\ &+ \frac{\textcolor{green}{{p}_{\rm5}(V_{\rm M}, d_{\rm r}, d_{\rm l})}}{\left(\textcolor{green}{{p}_{\rm6}(V_{\rm M},d_r,d_l)}+\textcolor{green}{{p}_{\rm7}(V_{\rm M},d_r,d_l)}\cdot\left(\textcolor{green}{{p}_{\rm8}(V_{\rm M},d_r,d_l)}e^{\ln\left(\frac{N_{\rm d}}{N_{\rm d,L}}\right)-\textcolor{green}{{p}_{\rm9}(V_{\rm M},d_r,d_l)}}\right)^{-\textcolor{green}{{p}_{\rm10}(V_{\rm M},d_r,d_l)}}\right)^{1/\textcolor{green}{{p}_{\rm11}(V_{\rm M},d_r,d_l)}}}
  \end{aligned}
  ```
#### $p_i$ functions
##### Applied Voltage V<sub>M</sub>>0
$p_1(V_{\rm M}, d_{\rm r}, d_{\rm l}) = 0$

$p_2(V_{\rm M}, d_{\rm r}, d_{\rm l}) = 0$


$p_3(V_{\rm M}, d_{\rm r}, d_{\rm l}) = 0$

$p_4(V_{\rm M}, d_{\rm r}, d_{\rm l}) = 1$


$p_5(V_{\rm M}, d_{\rm r}, d_{\rm l}) = \color{blue}p_{5,0|f}(d_{\rm r}, d_{\rm l})\color{black} - p_{5,1|f}(d_{\rm r}, d_{\rm l}) * e^{-p_{5,2|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}}$


```math
\begin{aligned}[b]
p_6(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{6,0|f}(d_{\rm r}, d_{\rm l}) + p_{6,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}
\end{aligned}
```


```math
\begin{aligned}[b]
p_7(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{7,0|f}(d_{\rm r}, d_{\rm l}) + p_{7,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M} + p_{7,2|f}(d_{\rm r}, d_{\rm l}) * e^{-p_{7,3|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}}
\end{aligned}
```


```math
\begin{aligned}[b]
p_8(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{8,0|f}(d_{\rm r}, d_{\rm l}) + p_{8,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}
\end{aligned}
```


```math
\begin{aligned}[b]
p_9(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & 0
\end{aligned}
```


```math
\begin{aligned}[b]
p_{10}(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{10,0|f}(d_{\rm r}, d_{\rm l}) + p_{10,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M} + p_{10,2|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}^2
\end{aligned}
```


```math
\begin{aligned}[b]
p_{11}(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{11,0|f}(d_{\rm r}, d_{\rm l}) + p_{11,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M} + p_{11,2|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}^2
\end{aligned}
```

##### Applied Voltage V<sub>M</sub><0
```math
\begin{aligned}[b]
p_1(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{1,0|f}(d_{\rm r}, d_{\rm l})\frac{p_{1,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M} + p_{1,2|f}(d_{\rm r}, d_{\rm l}) * V_{\rm M}^2}{1 + p_{1,3|f}(d_{\rm r}, d_{\rm l})*V_{\rm M} + p_{1,4|f}(d_{\rm r}, d_{\rm l}) * V_{\rm M}^2}\\
p_2(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{2,0|f}(d_{\rm r}, d_{\rm l})\\
p_3(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{3,0|f}(d_{\rm r}, d_{\rm l}) + p_{3,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}\\
p_4(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{4,0|f}(d_{\rm r}, d_{\rm l}) - p_{4,1|f}(d_{\rm r}, d_{\rm l}) * e^{-p_{4,2|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}}\\
p_5(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{5,0|f}(d_{\rm r}, d_{\rm l}) + p_{5,1|f}(d_{\rm r}, d_{\rm l})*V_{\rm M} + p_{5,2|f}(d_{\rm r}, d_{\rm l})*V_{\rm M}^2\\
p_6(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & 1\\
p_7(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{7,0|f}(d_{\rm r}, d_{\rm l})\\
p_8(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & 1
\end{aligned}
```
```math
\begin{aligned}[b]
p_9(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{9,0|f}(d_{\rm r}, d_{\rm l}) + \frac{p_{9,1|f}(d_{\rm r}, d_{\rm l}) - p_{9,0|f}(d_{\rm r}, d_{\rm l})}{1 + e^{\frac{V_{\rm M}-p_{9,2|f}(d_{\rm r}, d_{\rm l})}{p_{9,3|f}(d_{\rm r}, d_{\rm l})}}}\\
\end{aligned}
```

```math
\begin{aligned}[b]
p_{10}(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ & p_{10,0|f}(d_{\rm r}, d_{\rm l}) + \frac{p_{10,1|f}(d_{\rm r}, d_{\rm l}) - p_{10,0|f}(d_{\rm r}, d_{\rm l})}{1 + e^{\frac{V_{\rm M}-p_{10,2|f}(d_{\rm r}, d_{\rm l})}{p_{10,3|f}(d_{\rm r}, d_{\rm l})}}}\\
p_{11}(V_{\rm M}, d_{\rm r}, d_{\rm l}) ~=~ &p_{11,0|f}(d_{\rm r}, d_{\rm l}) + \frac{p_{11,1|f}(d_{\rm r}, d_{\rm l}) - p_{11,0|f}(d_{\rm r}, d_{\rm l})}{1 + e^{\frac{V_{\rm M}-p_{11,2|f}(d_{\rm r}, d_{\rm l})}{p_{11,3|f}(d_{\rm r}, d_{\rm l})}}}
\end{aligned}
```

## Fitting parameter values
### Applied Voltage V<sub>M</sub>>0
<table>
    <thead>
        <tr>
            <th>Parameter</th>
            <th>Vm</th>
            <th>d_r</th>
            <th>d_l</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>p<sub>1</sub></td>
            <td>0</td>
            <td></td>
            <td></td>
        </tr>
        <tr>
            <td>p<sub>2</sub></td>
            <td>0</td>
            <td></td>
            <td></td>
        </tr>
        <tr>
            <td>p<sub>3</sub></td>
            <td>0</td>
            <td></td>
            <td></td>
        </tr>
        <tr>
            <td>p<sub>4</sub></td>
            <td>1</td>
            <td></td>
            <td></td>
        </tr>
        <tr>
            <td>p<sub>5</sub></td>
            <td>
                p<sub>5,1</sub>=1.3769e-03 <br>
                p<sub>5,2</sub>=8.1819e-02
            </td>
            <td>
                D<sub>p<sub>5,1,r</sub></sub>=2.3087e-04
            </td>
            <td>
                D<sub>p<sub>5,1,l</sub></sub>=1.3293e-07
            </td>
        </tr>
        <tr>
            <td>p<sub>6</sub></td>
            <td>
                p<sub>6,0</sub>=1.9687e-01 <br>
                p<sub>6,1</sub>=-2.1833e-02
            </td>
            <td>
                D<sub>p<sub>6,0,r</sub></sub>=2.6129e-02
            </td>
            <td>
            </td>
        </tr>
        <tr>
            <td>p<sub>7</sub></td>
            <td>
                p<sub>7,0</sub>=-9.7606e+01 <br>
                p<sub>7,1</sub>=7.8250e+00 <br>
                p<sub>7,2</sub>=9.9296e+01 <br>
                p<sub>7,3</sub>=7.1092e-02
            </td>
            <td>
                D<sub>p<sub>7,0,r</sub></sub>=-7.4338e-01 <br/>
                D<sub>p<sub>7,0,r<sup>2</sup></sub></sub>=1.1713e-01 <br/>
                D<sub>p<sub>7,2,r</sub></sub>=6.9547e-01 <br/>
                D<sub>p<sub>7,2,r<sup>2</sup></sub></sub>=-1.3724e-01 <br/>
                D<sub>p<sub>7,2,l,r</sub></sub>=-9.0456e-03 <br/>
                D<sub>p<sub>7,2,l,r<sup>2</sup></sub></sub>=-1.2221e-03
            </td>
            <td>
                D<sub>p<sub>7,0,l</sub></sub>=2.4377e+00 <br/>
                D<sub>p<sub>7,2,l</sub></sub>=-2.3728e+00
            </td>
        </tr>
        <tr>
            <td>p<sub>8</sub></td>
            <td>
                p<sub>8,0</sub>=1.1713e-01 <br>
                p<sub>8,1</sub>=8.1370e-02
            </td>
            <td></td>
            <td>
                D<sub>p<sub>8,0,l</sub></sub>=-3.8320e-03
            </td>
        </tr>
        <tr>
            <td>p<sub>9</sub></td>
            <td>0</td>
            <td></td>
            <td></td>
        </tr>
        <tr>
            <td>p<sub>10</sub></td>
            <td>
                p<sub>10,0</sub>=9.7733e-01 <br>
                p<sub>10,1</sub>=3.5214e-02 <br>
                p<sub>10,2</sub>=1.2856e-02
            </td>
            <td></td>
            <td>
                D<sub>p<sub>10,0,l</sub></sub>=5.9623e-05
            </td>
        </tr>
        <tr>
            <td>p<sub>11</sub></td>
            <td>
                p<sub>11,0</sub>=9.4207e-01 <br>
                p<sub>11,1</sub>=3.8953e-02 <br>
                p<sub>11,2</sub>=2.3436e-02
            </td>
            <td></td>
            <td>
                D<sub>p<sub>11,0,l</sub></sub>=-6.7239e-04
            </td>
        </tr>
    </tbody>
</table>

### Applied Voltage V<sub>M</sub><0
<table>
    <thead>
        <tr>
            <th>Parameter</th>
            <th>Vm</th>
            <th>d_r</th>
            <th>d_l</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>p<sub>1</sub></td>
            <td>
                p<sub>1,0</sub>=1.1830e+00 <br>
                p<sub>1,1</sub>=-2.7034e-03 <br>
                p<sub>1,2</sub>=-4.5379e-06 <br>
                p<sub>1,3</sub>=9.9115e-01 <br>
                p<sub>1,4</sub>=4.4093e-01
            </td>
            <td>
                D<sub>p<sub>1,0,r</sub></sub>=-6.2246e-02 <br>
                D<sub>p<sub>1,1,r</sub></sub>=-1.8077e-04 <br>
                D<sub>p<sub>1,2,r</sub></sub>=1.7313e-04 <br>
                D<sub>p<sub>1,3,r</sub></sub>=5.7155e-03 <br>
                D<sub>p<sub>1,4,r</sub></sub>=-9.7198e-04
            </td>
            <td>
                D<sub>p<sub>1,0,l</sub></sub>=1.1419e-01 <br>
                D<sub>p<sub>1,1,l</sub></sub>=-2.0831e-04 <br>
                D<sub>p<sub>1,2,l</sub></sub>=-8.9677e-05 <br>
                D<sub>p<sub>1,3,l</sub></sub>=-2.3237e-02 <br>
                D<sub>p<sub>1,4,l</sub></sub>=-1.8507e-03
            </td>
        </tr>
        <tr>
            <td>p<sub>2</sub></td>
            <td>-2.5955e+03</td>
            <td></td>
            <td></td>
        </tr>
        <tr>
            <td>p<sub>3</sub></td>
            <td>
                p<sub>3,0</sub>=6.8845e+00 <br>
                p<sub>3,1</sub>=-5.8995e-01
            </td>
            <td>
                D<sub>p<sub>3,0,r</sub></sub>=1.2536e-01 <br>
                D<sub>p<sub>3,1,r</sub></sub>=6.5498e-02
            </td>
            <td>
                D<sub>p<sub>3,0,l</sub></sub>=2.5983e-01 <br>
                D<sub>p<sub>3,1,l</sub></sub>=8.5666e-02
            </td>
        </tr>
        <tr>
            <td>p<sub>4</sub></td>
            <td>
                p<sub>4,0</sub>=2.5890e+03 <br>
                p<sub>4,1</sub>=-2.9537e+00 <br>
                p<sub>4,2</sub>=-5.4031e-01
            </td>
            <td>
                D<sub>p<sub>4,1,r</sub></sub>=8.2522e-02
            </td>
            <td>
                D<sub>p<sub>4,1,l</sub></sub>=-7.2255e-02
            </td>
        </tr>
        <tr>
            <td>p<sub>5</sub></td>
            <td>
                p<sub>5,0</sub>=0 <br>
                p<sub>5,1</sub>=6.4705e-04 <br>
                p<sub>5,2</sub>=5.1529e-05
            </td>
            <td>
                D<sub>p<sub>5,1,r</sub></sub>=1.5169e-05 <br>
                D<sub>p<sub>5,2,r</sub></sub>=6.7042e-07 <br>
                D<sub>p<sub>5,2,r<sup>2</sup></sub></sub>=1.0756e-06
            </td>
            <td>
                D<sub>p<sub>5,1,l</sub></sub>=1.3260e-06
            </td>
        </tr>
        <tr>
            <td>p<sub>6</sub></td>
            <td>
                1
            </td>
            <td>
            </td>
            <td>
            </td>
        </tr>
        <tr>
            <td>p<sub>7</sub></td>
            <td>
                1.1708e-01
            </td>
            <td>
                D<sub>p<sub>7,0,r</sub></sub>=4.8662e-04
            </td>
            <td>
                D<sub>p<sub>7,0,l</sub></sub>=3.7351e-03
            </td>
        </tr>
        <tr>
            <td>p<sub>8</sub></td>
            <td>
                1
            </td>
            <td></td>
            <td>
            </td>
        </tr>
        <tr>
            <td>p<sub>9</sub></td>
            <td>
                p<sub>9,0</sub>=3.9052e+00 <br>
                p<sub>9,1</sub>=9.6130e+00 <br>
                p<sub>9,2</sub>=-4.5637e-01 <br>
                p<sub>9,3</sub>=1.4310e+00
            </td>
            <td>
                D<sub>p<sub>9,0,r</sub></sub>=-5.4723e-01 <br/>
                D<sub>p<sub>9,3,r</sub></sub>=3.6000e-01
            </td>
            <td>
                D<sub>p<sub>9,0,l</sub></sub>=3.6802e-02
            </td>
        </tr>
        <tr>
            <td>p<sub>10</sub></td>
            <td>
                p<sub>10,0</sub>=4.6925e-01 <br>
                p<sub>10,1</sub>=3.4731e+00 <br>
                p<sub>10,2</sub>=-1.1871e+00 <br>
                p<sub>10,3</sub>=5.6947e-01
            </td>
            <td>
                D<sub>p<sub>10,1,r</sub></sub>=1.1444e-02
            </td>
            <td>
            </td>
        </tr>
        <tr>
            <td>p<sub>11</sub></td>
            <td>
                p<sub>11,0</sub>=1.0667e+01 <br>
                p<sub>11,1</sub>=1.2812e-01 <br>
                p<sub>11,2</sub>=7.4414e-01 <br>
                p<sub>11,3</sub>=4.2381e-01 
            </td>
            <td>
                D<sub>p<sub>11,0,r</sub></sub>=3.6290e-01
            </td>
            <td>
            </td>
        </tr>
    </tbody>
</table>
