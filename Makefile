CC		= clang++
CPPFLAGS	= -Wall -Wextra -std=c++1y
DBGTARGET  	= ${TARGET}_dbg
DBGOBJS		= $(patsubst ${SRCDIR}%,${OBJDIR}%,${SRCS:.cc=.dbgo})
DBGFLAGS	= -O0 -ggdb3
OBJDIR		= objects
OBJS		= $(patsubst ${SRCDIR}%,${OBJDIR}%,${SRCS:.cc=.o})
OPTFLAGS	= -O3 -DNDEBUG
SRCDIR		= src
SRCS		= $(wildcard ${SRCDIR}/*.cc)
TARGET		= heat

all: ${TARGET}

${OBJDIR}/:
	mkdir -p ${OBJDIR}

${OBJDIR}/%.dbgo : ${SRCDIR}/%.cc | ${OBJDIR}/
	${CC} ${CPPFLAGS} ${DBGFLAGS} -c $< -o $@

${OBJDIR}/%.o : ${SRCDIR}/%.cc | ${OBJDIR}/
	${CC} ${CPPFLAGS} ${OPTFLAGS} -c $< -o $@

${TARGET}: ${OBJS}
	${CC} ${OBJS} -o ${TARGET}

${DBGTARGET}: ${DBGOBJS}
	${CC} ${DBGOBJS} ${DBGFLAGS} -o ${DBGTARGET}

clean:
	rm -rf ${TARGET} ${DBGTARGET} ${OBJDIR} core

debug: ${DBGTARGET}
	gdb ${DBGTARGET} -x debug_script

debug_auto: ${DBGTARGET}
	gdb ${DBGTARGET} -x debug_script_auto

run: ${DBGTARGET}
	./${DBGTARGET}