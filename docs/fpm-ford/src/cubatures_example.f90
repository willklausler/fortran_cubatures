program cubatures_example
!! Demonstrate cubatures

  use cubatures

  implicit none

  type(Cubature) :: scheme

  call scheme%set("LIN",[2])

  call scheme%show()

  call scheme%set("TRI",[2])

  call scheme%show()

  call scheme%set("QUA",[2])

  call scheme%show()

  call scheme%set("TET",[2])

  call scheme%show()

  call scheme%set("HEX",[2])

  call scheme%show()

  call scheme%set("WEJ",[1])

  call scheme%show()

end program cubatures_example