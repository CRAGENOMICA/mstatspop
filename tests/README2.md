# mstatspop

## Testing
```
## npm is required
sudo npm install -g bats

npm install  https://github.com/ztombol/bats-assert
npm install  https://github.com/ztombol/bats-support
npm install  https://github.com/ztombol/bats-file
```


# To Run the Test
```
bats ./tests/mstatspop.tests.bats
```

# Sanityzing

```
BUILD_TYPE=Debug CFLAGS="-fsanitize=address -g" sh build.sh
#BUILD_TYPE=Release CFLAGS="-O3 -fsanitize=address -g" sh build.sh

bats ./tests/mstatspop.tests.bats
```
