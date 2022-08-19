# ECM-Rust
楕円曲線法を使った素因数分解プログラムです．

## 引数
```
factorize a number using elliptic curve method

USAGE:
    ecm_rust [OPTIONS] --n <n> --b1 <b1> --b2 <b2> --d <d> --k <k> --thread_num <thread_num>

OPTIONS:
    -b, --b1 <b1>                    the bound of primes in the step 1
    -B, --b2 <b2>                    the bound of primes in the step 2
    -d, --d <d>                      the chunk size of the step 2
    -h, --help                       Print help information
    -k, --k <k>                      the degree of power used in the step 2
    -n, --n <n>                      the number to factorize
    -s, --seed <seed>                the seed to use
    -t, --thread_num <thread_num>    the number of threads to use
    -V, --version                    Print version information
```


## 使用例

```
$ cargo build --release
$ ./target/release/ecm_rust -n 231584178474632390847141970017375815706539969331281128078915168015826259279871 \
                            -b 300000 -B 10000000 -d 4096 -t 12 -k 12 -s 1234
seed : 1234

factorizing : 231584178474632390847141970017375815706539969331281128078915168015826259279871
factorize attempt count : 12 (total : 12)  time : 2.887653412s
factors found :
231584178474632390847141970017375815706539969331281128078915168015826259279871
         = 535006138814359 (15 digits, prime)
         x 432862656469423142931042426214547535783388063929571229938474969 (63 digits, composite)
sigma : [1358743495426659592, 3702083616141384991, 5057458246547433278, 347741063316059062, 6036757013344071376, 8256042427999371329, 4278107792852053917]

factorizing : 432862656469423142931042426214547535783388063929571229938474969
factorize attempt count : 12 (total : 24)  time : 2.505255186s
factorize attempt count : 24 (total : 36)  time : 2.452675097s
factors found :
432862656469423142931042426214547535783388063929571229938474969
         = 1155685395246619182673033 (25 digits, prime)
         x 374550598501810936581776630096313181393 (39 digits, prime)
sigma : [967271087794501515]

result :
231584178474632390847141970017375815706539969331281128078915168015826259279871
         = 535006138814359 (15 digits, prime)
         x 1155685395246619182673033 (25 digits, prime)
         x 374550598501810936581776630096313181393 (39 digits, prime)
process time: 9.012141411s
```

## パラメータの決定
パラメータは見つけたい素因子の桁数によっておおよそ決まります．  

b1, b2 はそれぞれ step1, step2 で使用する素数の上界です．  
b1, b2 ともに値を増やすと一回の施行で成功する確率は増えますが，処理時間は毎回おおよそ O(b1 + √b2) となるので増やしすぎると逆効果です．  
経験的には b1 は見つけたい素因子の桁数を 30, 40, 50 桁とするとそれぞれ 3e5, 3e6 5e7 程度，b2 は b1 の 1000 ~ 100000 倍にすると良いです．  
この値については GMP-ECM の資料が参考になります． ( https://github.com/sethtroisi/gmp-ecm より引用)

```
       digits D  optimal B1   default B2           expected curves
                                                       N(B1,B2,D)
                                              -power 1         default poly
          20       11e3         1.9e6             74               74 [x^1]
          25        5e4         1.3e7            221              214 [x^2]
          30       25e4         1.3e8            453              430 [D(3)]
          35        1e6         1.0e9            984              904 [D(6)]
          40        3e6         5.7e9           2541             2350 [D(6)]
          45       11e6        3.5e10           4949             4480 [D(12)]
          50       43e6        2.4e11           8266             7553 [D(12)]
          55       11e7        7.8e11          20158            17769 [D(30)]
          60       26e7        3.2e12          47173            42017 [D(30)]
          65       85e7        1.6e13          77666            69408 [D(30)]
```

d は step2 で一度に計算されるチャンクの大きさです．  
基本的には　√b2 程度の二冪にすると最適ですが，あまりに大きいとメモリを食います．

k は step2 で使用する 2 と 3 の冪の次数です．  
k の約数が多いほど効率的になりますが，あまりに k が大き過ぎると step2 に k に比例する時間がかかるようになるので動作が遅くなります．  
素因数の大きさに応じて 12, 60, 120, 180, 240 ぐらいの値を使い分けると良いと思います．
