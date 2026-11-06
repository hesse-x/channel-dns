! ---&---1---------2---------3---------4---------5---------6---------7--
! Configuration subset: flat sections, scalar values and numeric arrays.
MODULE dns_toml
  USE dns_kinds, ONLY: dp, INT64, fail
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: IOSTAT_END, IOSTAT_EOR
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: toml_config, toml_table, table, has_key, check_keys
  PUBLIC :: read_integer, read_real, read_vector, read_string
  TYPE :: ENTRY
    CHARACTER(:), ALLOCATABLE :: key, VALUE
  END TYPE
  TYPE :: toml_table
    CHARACTER(:), ALLOCATABLE :: name
    TYPE(ENTRY), ALLOCATABLE :: items(:)
    TYPE(toml_table), ALLOCATABLE :: children(:)
  END TYPE
  TYPE :: toml_config
    TYPE(toml_table), ALLOCATABLE :: handle
  CONTAINS
    PROCEDURE :: load
    PROCEDURE :: CLOSE
  END TYPE
CONTAINS

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  FUNCTION lower(text) RESULT(VALUE)
    CHARACTER(*), INTENT(in) :: text
    CHARACTER(LEN(text)) :: VALUE
    INTEGER :: i, code

    VALUE = text
    DO i = 1, LEN(text)
      code = IACHAR(text(i:i))
      IF (code >= 65 .AND. code <= 90) VALUE(i:i) = ACHAR(code+32)
    END DO
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  LOGICAL FUNCTION valid_name(name)
    CHARACTER(*), INTENT(in) :: name

    valid_name = LEN(name) > 0
    IF (VERIFY(name, 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-') /= 0) valid_name = .FALSE.
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE syntax_error(path, line, message)
    CHARACTER(*), INTENT(in) :: path, message
    INTEGER, INTENT(in) :: line
    CHARACTER(24) :: number

    WRITE (number, '(i0)') line
    CALL fail(TRIM(path)//':'//TRIM(number)//': '//message)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE load(self, path)
    CLASS(toml_config), TARGET, INTENT(inout) :: self
    CHARACTER(*), INTENT(in) :: path
    TYPE(toml_table), POINTER :: current
    CHARACTER(4097) :: buffer
    CHARACTER(:), ALLOCATABLE :: text, name, VALUE
    CHARACTER :: quote, ch
    INTEGER :: unit, ios, line, i, cut, eq, n, count

    CALL self%CLOSE()
    ALLOCATE (self%handle)
    self%handle%name = 'root'
    ALLOCATE (self%handle%items(0), self%handle%children(0))
    NULLIFY (current)
    OPEN (newunit=unit, file=path, status='old', action='read', iostat=ios)
    IF (ios /= 0) CALL fail('cannot open '//path)
    line = 0
    DO
      READ (unit, '(a)', advance='no', size=count, iostat=ios) buffer
      IF (ios == IOSTAT_END) EXIT
      line = line+1
      IF (count > 4096) CALL syntax_error(path, line, 'line exceeds 4096 characters')
      IF (ios /= IOSTAT_EOR) CALL syntax_error(path, line, 'cannot read configuration')
      quote = ' '
      cut = LEN_TRIM(buffer)+1
      DO i = 1, LEN_TRIM(buffer)
        ch = buffer(i:i)
        IF (quote == ' ' .AND. ch == ACHAR(9)) buffer(i:i) = ' '
        IF (quote /= ' ') THEN
          IF (ch == quote) quote = ' '
        ELSE IF (ch == '"' .OR. ch == "'") THEN
          quote = ch
        ELSE IF (ch == '#') THEN
          cut = i
          EXIT
        END IF
      END DO
      IF (quote /= ' ') CALL syntax_error(path, line, 'unterminated string')
      text = TRIM(ADJUSTL(buffer(:cut-1)))
      IF (LEN(text) == 0) CYCLE
      IF (text(1:1) == '[') THEN
        n = LEN(text)
        IF (text(n:n) /= ']') CALL syntax_error(path, line, 'invalid section')
        name = lower(TRIM(ADJUSTL(text(2:n-1))))
        IF (.NOT. valid_name(name)) CALL syntax_error(path, line, 'simple section name required')
        IF (has_key(self%handle, name)) CALL syntax_error(path, line, 'duplicate section: '//name)
        NULLIFY (current)
        CALL add_section(self%handle, name)
        n = SIZE(self%handle%children)
        current => self%handle%children(n)
      ELSE
        IF (.NOT. ASSOCIATED(current)) CALL syntax_error(path, line, 'key must belong to a section')
        eq = INDEX(text, '=')
        IF (eq < 2) CALL syntax_error(path, line, 'key = value required')
        name = lower(TRIM(text(:eq-1)))
        IF (.NOT. valid_name(name)) CALL syntax_error(path, line, 'simple key name required')
        IF (has_key(current, name)) CALL syntax_error(path, line, 'duplicate key: '//name)
        VALUE = TRIM(ADJUSTL(text(eq+1:)))
        IF (.NOT. valid_value(VALUE)) CALL syntax_error(path, line, 'unsupported or invalid value for '//name)
        CALL add_entry(current, name, VALUE)
      END IF
    END DO
    CLOSE (unit)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE add_section(root, name)
    TYPE(toml_table), INTENT(inout) :: root
    CHARACTER(*), INTENT(in) :: name
    TYPE(toml_table), ALLOCATABLE :: grown(:)
    INTEGER :: i, n

    n = SIZE(root%children)
    ALLOCATE (grown(n+1))
    DO i = 1, n
      CALL MOVE_ALLOC(root%children(i)%name, grown(i)%name)
      CALL MOVE_ALLOC(root%children(i)%items, grown(i)%items)
      CALL MOVE_ALLOC(root%children(i)%children, grown(i)%children)
    END DO
    grown(n+1)%name = name
    ALLOCATE (grown(n+1)%items(0), grown(n+1)%children(0))
    CALL MOVE_ALLOC(grown, root%children)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE add_entry(section, name, VALUE)
    TYPE(toml_table), INTENT(inout) :: section
    CHARACTER(*), INTENT(in) :: name, VALUE
    TYPE(ENTRY), ALLOCATABLE :: grown(:)
    INTEGER :: i, n

    n = SIZE(section%items)
    ALLOCATE (grown(n+1))
    DO i = 1, n
      CALL MOVE_ALLOC(section%items(i)%key, grown(i)%key)
      CALL MOVE_ALLOC(section%items(i)%VALUE, grown(i)%VALUE)
    END DO
    grown(n+1)%key = name
    grown(n+1)%VALUE = VALUE
    CALL MOVE_ALLOC(grown, section%items)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE CLOSE (self)
    CLASS(toml_config), INTENT(inout) :: self

    IF (ALLOCATED(self%handle)) DEALLOCATE (self%handle)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  FUNCTION table(handle, key) RESULT(child)
    TYPE(toml_table), TARGET, INTENT(inout) :: handle
    CHARACTER(*), INTENT(in) :: key
    TYPE(toml_table), POINTER :: child
    INTEGER :: i

    NULLIFY (child)
    DO i = 1, SIZE(handle%children)
      IF (handle%children(i)%name /= lower(key)) CYCLE
      child => handle%children(i)
      RETURN
    END DO
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  LOGICAL FUNCTION has_key(handle, key)
    TYPE(toml_table), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: key
    INTEGER :: i

    has_key = .TRUE.
    DO i = 1, SIZE(handle%children)
      IF (handle%children(i)%name == lower(key)) RETURN
    END DO
    DO i = 1, SIZE(handle%items)
      IF (handle%items(i)%key == lower(key)) RETURN
    END DO
    has_key = .FALSE.
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE check_keys(handle, allowed, path)
    TYPE(toml_table), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: allowed, path
    CHARACTER(:), ALLOCATABLE :: key
    INTEGER :: i

    DO i = 1, SIZE(handle%children)+SIZE(handle%items)
      IF (i <= SIZE(handle%children)) THEN
        key = handle%children(i)%name
      ELSE
        key = handle%items(i-SIZE(handle%children))%key
      END IF
      IF (INDEX('|'//allowed//'|', '|'//key//'|') == 0) CALL fail('unknown TOML name: '//path//'.'//key)
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  FUNCTION raw_value(handle, key) RESULT(VALUE)
    TYPE(toml_table), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: key
    CHARACTER(:), ALLOCATABLE :: VALUE
    INTEGER :: i

    DO i = 1, SIZE(handle%items)
      IF (handle%items(i)%key /= lower(key)) CYCLE
      VALUE = handle%items(i)%VALUE
      RETURN
    END DO
    CALL fail('missing TOML key: '//handle%name//'.'//key)
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
! Return 1 for decimal integers, 2 for reals and 0 for invalid tokens.
  INTEGER FUNCTION number_kind(text) RESULT(kind)
    CHARACTER(*), INTENT(in) :: text
    INTEGER :: i, first, n

    kind = 0
    n = LEN(text)
    IF (n == 0) RETURN
    i = 1
    IF (text(i:i) == '+' .OR. text(i:i) == '-') i = i+1
    first = i
    DO WHILE (i <= n)
      IF (VERIFY(text(i:i), '0123456789') /= 0) EXIT
      i = i+1
    END DO
    IF (i == first) RETURN
    IF (i-first > 1) THEN
      IF (text(first:first) == '0') RETURN
    END IF
    kind = 1
    IF (i <= n) THEN
      IF (text(i:i) == '.') THEN
        kind = 2
        i = i+1
        first = i
        DO WHILE (i <= n)
          IF (VERIFY(text(i:i), '0123456789') /= 0) EXIT
          i = i+1
        END DO
        IF (i == first) kind = 0
      END IF
    END IF
    IF (i <= n) THEN
      IF (text(i:i) == 'e' .OR. text(i:i) == 'E') THEN
        IF (kind == 0) RETURN
        kind = 2
        i = i+1
        IF (i <= n) THEN
          IF (text(i:i) == '+' .OR. text(i:i) == '-') i = i+1
        END IF
        first = i
        DO WHILE (i <= n)
          IF (VERIFY(text(i:i), '0123456789') /= 0) EXIT
          i = i+1
        END DO
        IF (i == first) kind = 0
      END IF
    END IF
    IF (i <= n) kind = 0
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  LOGICAL FUNCTION valid_value(VALUE) RESULT(valid)
    CHARACTER(*), INTENT(in) :: VALUE
    INTEGER :: n, first, pos, next

    valid = .FALSE.
    n = LEN(VALUE)
    IF (n == 0) RETURN
    IF (VALUE(1:1) == '"' .OR. VALUE(1:1) == "'") THEN
      IF (n < 2) RETURN
      IF (VALUE(n:n) /= VALUE(1:1)) RETURN
      IF (INDEX(VALUE(2:n-1), VALUE(1:1)) /= 0) RETURN
      IF (INDEX(VALUE, ACHAR(92)) /= 0) RETURN
      valid = .TRUE.
    ELSE IF (VALUE(1:1) == '[') THEN
      IF (VALUE(n:n) /= ']') RETURN
      first = 2
      DO
        pos = INDEX(VALUE(first:n-1), ',')
        next = n
        IF (pos > 0) next = first+pos-1
        IF (number_kind(TRIM(ADJUSTL(VALUE(first:next-1)))) == 0) RETURN
        IF (pos == 0) EXIT
        first = next+1
      END DO
      valid = .TRUE.
    ELSE
      valid = number_kind(VALUE) /= 0
    END IF
  END FUNCTION

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE read_integer(handle, key, VALUE)
    TYPE(toml_table), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: key
    INTEGER, INTENT(out) :: VALUE
    CHARACTER(:), ALLOCATABLE :: text
    INTEGER(INT64) :: wide
    INTEGER :: ios

    text = raw_value(handle, key)
    IF (number_kind(text) /= 1) CALL fail('integer required: '//key)
    READ (text, *, iostat=ios) wide
    IF (ios /= 0) CALL fail('integer overflow: '//key)
    IF (wide > INT(HUGE(VALUE), INT64) .OR. wide < (-INT(HUGE(VALUE), INT64)-1)) CALL fail('integer overflow: '//key)
    VALUE = INT(wide)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE read_real(handle, key, VALUE)
    TYPE(toml_table), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: key
    REAL(dp), INTENT(out) :: VALUE
    CHARACTER(:), ALLOCATABLE :: text
    INTEGER :: ios

    text = raw_value(handle, key)
    IF (number_kind(text) == 0) CALL fail('number required: '//key)
    READ (text, *, iostat=ios) VALUE
    IF (ios /= 0) CALL fail('number out of range: '//key)
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE read_vector(handle, key, values)
    TYPE(toml_table), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: key
    REAL(dp), INTENT(out) :: values(3)
    CHARACTER(:), ALLOCATABLE :: text, item
    INTEGER :: i, first, pos, next, n, ios

    text = raw_value(handle, key)
    n = LEN(text)
    IF (text(1:1) /= '[') CALL fail('numeric array required: '//key)
    first = 2
    DO i = 1, 3
      pos = INDEX(text(first:n-1), ',')
      next = n
      IF (pos > 0) next = first+pos-1
      IF ((i < 3 .AND. pos == 0) .OR. (i == 3 .AND. pos /= 0)) CALL fail('three-element array required: '//key)
      item = TRIM(ADJUSTL(text(first:next-1)))
      READ (item, *, iostat=ios) values(i)
      IF (ios /= 0) CALL fail('number out of range: '//key)
      first = next+1
    END DO
  END SUBROUTINE

! ---&---1---------2---------3---------4---------5---------6---------7--
!
  SUBROUTINE read_string(handle, key, VALUE)
    TYPE(toml_table), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: key
    CHARACTER(:), ALLOCATABLE, INTENT(out) :: VALUE
    CHARACTER(:), ALLOCATABLE :: text
    INTEGER :: n

    text = raw_value(handle, key)
    n = LEN(text)
    IF (text(1:1) /= '"' .AND. text(1:1) /= "'") CALL fail('string required: '//key)
    VALUE = text(2:n-1)
  END SUBROUTINE
END MODULE
! ---&---1---------2---------3---------4---------5---------6---------7--
