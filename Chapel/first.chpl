var sum : int = 0;
forall i in 1..1000000 with (ref sum) {
    sum += i; //DANGER: race condition!
}
writeln(sum);
